/***********************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2026, Meng Wang                                        *
***********************************************************************/

#include <cuda_runtime.h>

#include "lucia_cuda_buffer.cuh"
#include "lucia_sigma_cuda_blocks.cuh"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>

namespace {

constexpr std::size_t tile_size = 16;
constexpr std::size_t scale_block_size = 256;
constexpr std::size_t reduction_target_per_split = 1024;
constexpr std::size_t max_split_count = 16;

bool checked_add(std::size_t left, std::size_t right, std::size_t *result) noexcept {
  if (left > (std::numeric_limits<std::size_t>::max)() - right) {
    return false;
  }
  *result = left + right;
  return true;
}
bool checked_triangular(std::size_t n, std::size_t *result) noexcept {
  std::size_t left = 0;
  std::size_t right = 0;
  if ((n & 1U) == 0U) {
    left = n / 2;
    if (!checked_add(n, 1, &right)) {
      return false;
    }
  } else {
    left = n;
    right = n / 2 + 1;
  }
  return lucia_cuda::checked_mul(left, right, result);
}

bool append_section(std::size_t bytes, std::size_t *offset, std::size_t *total) noexcept {
  *offset = *total;
  return checked_add(*total, bytes, total);
}

bool checked_ceil_div(std::size_t numerator, std::size_t denominator, std::size_t *result) noexcept {
  if (denominator == 0) {
    return false;
  }
  *result = numerator / denominator;
  if (numerator % denominator != 0 && !checked_add(*result, 1, result)) {
    return false;
  }
  return true;
}

bool valid_extent(std::size_t offset, std::size_t extent, std::size_t limit) noexcept {
  return offset > 0 && extent > 0 && offset <= limit && extent <= limit - (offset - 1);
}

bool fill_pair_columns(std::int64_t flag, std::size_t first_extent, std::size_t second_extent, std::size_t pair_count,
                       std::int64_t *columns) noexcept {
  if (columns == nullptr || pair_count == 0 || pair_count > static_cast<std::size_t>((std::numeric_limits<std::int64_t>::max)())) {
    return false;
  }

  if (flag == 0) {
    std::size_t expected = 0;
    if (!lucia_cuda::checked_mul(first_extent, second_extent, &expected) || expected != pair_count) {
      return false;
    }
    for (std::size_t p = 0; p < pair_count; ++p) {
      columns[p] = static_cast<std::int64_t>(p);
    }
    return true;
  }

  if (flag != 1 || first_extent != second_extent) {
    return false;
  }
  std::size_t expected = 0;
  if (!checked_triangular(first_extent, &expected) || expected != pair_count) {
    return false;
  }

  std::size_t first0 = 0;
  std::size_t second0 = 0;
  for (std::size_t p = 0; p < pair_count; ++p) {
    if (first0 == second0) {
      columns[p] = -1;
    } else {
      std::size_t product = 0;
      std::size_t full_column = 0;
      if (!lucia_cuda::checked_mul(second0, first_extent, &product) || !checked_add(product, first0, &full_column)
          || full_column > static_cast<std::size_t>((std::numeric_limits<std::int64_t>::max)())) {
        return false;
      }
      columns[p] = static_cast<std::int64_t>(full_column);
    }

    if (p < pair_count - 1) {
      if (second0 == first0) {
        if (!checked_add(first0, 1, &first0)) {
          return false;
        }
        second0 = 0;
      } else if (!checked_add(second0, 1, &second0)) {
        return false;
      }
    }
  }
  return true;
}

bool validate_map_columns(const std::int64_t *maps, std::int64_t flag, std::size_t first_extent, std::size_t second_extent,
                          std::size_t pair_count, std::size_t maxk, std::size_t used_k, std::int64_t max_source) noexcept {
  if (maps == nullptr || pair_count == 0 || pair_count > static_cast<std::size_t>((std::numeric_limits<std::int64_t>::max)())) {
    return false;
  }

  std::size_t expected = 0;
  if (flag == 0) {
    if (!lucia_cuda::checked_mul(first_extent, second_extent, &expected) || expected != pair_count) {
      return false;
    }
  } else if (flag != 1 || first_extent != second_extent || !checked_triangular(first_extent, &expected) || expected != pair_count) {
    return false;
  }

  const auto validate_column = [&](std::size_t column) noexcept {
    std::size_t map_base = 0;
    if (!lucia_cuda::checked_mul(column, maxk, &map_base)) {
      return false;
    }
    for (std::size_t q = 0; q < used_k; ++q) {
      std::size_t map_offset = 0;
      if (!checked_add(map_base, q, &map_offset) || maps[map_offset] < 0 || maps[map_offset] > max_source) {
        return false;
      }
    }
    return true;
  };

  if (flag == 0) {
    for (std::size_t column = 0; column < pair_count; ++column) {
      if (!validate_column(column)) {
        return false;
      }
    }
    return true;
  }

  std::size_t first0 = 0;
  std::size_t second0 = 0;
  for (std::size_t p = 0; p < pair_count; ++p) {
    if (first0 != second0) {
      std::size_t product = 0;
      std::size_t full_column = 0;
      if (!lucia_cuda::checked_mul(second0, first_extent, &product) || !checked_add(product, first0, &full_column)
          || full_column > static_cast<std::size_t>((std::numeric_limits<std::int64_t>::max)()) || !validate_column(full_column)) {
        return false;
      }
    }
    if (p < pair_count - 1) {
      if (second0 == first0) {
        if (!checked_add(first0, 1, &first0)) {
          return false;
        }
        second0 = 0;
      } else if (!checked_add(second0, 1, &second0)) {
        return false;
      }
    }
  }
  return true;
}

struct RouteSlot {
  lucia_cuda::DeviceBuffer<unsigned char> device;
  lucia_cuda::PinnedBuffer<unsigned char> host;
  cudaEvent_t host_ready = nullptr;
  bool pending = false;

  bool prepare() noexcept {
    if (!pending) {
      return true;
    }
    if (cudaEventSynchronize(host_ready) != cudaSuccess) {
      return false;
    }
    pending = false;
    return true;
  }

  bool ensure_event() noexcept {
    return host_ready != nullptr || cudaEventCreateWithFlags(&host_ready, cudaEventDisableTiming) == cudaSuccess;
  }

  bool upload(std::size_t bytes) noexcept {
    if (!device.reserve(bytes) || !ensure_event()
        || cudaMemcpyAsync(device.data, host.data, bytes, cudaMemcpyHostToDevice, 0) != cudaSuccess
        || cudaEventRecord(host_ready, 0) != cudaSuccess) {
      return false;
    }
    pending = true;
    return true;
  }

  bool release() noexcept {
    bool success = true;
    if (pending && cudaEventSynchronize(host_ready) != cudaSuccess) {
      success = false;
    }
    pending = false;
    if (host_ready != nullptr && cudaEventDestroy(host_ready) != cudaSuccess) {
      success = false;
    }
    host_ready = nullptr;
    if (!device.release())
      success = false;
    if (!host.release())
      success = false;
    return success;
  }
};

struct Workspace {
  lucia_cuda::DeviceBuffer<double> x;
  lucia_cuda::DeviceBuffer<double> sb;
  lucia_cuda::DeviceBuffer<double> cb;
  lucia_cuda::DeviceBuffer<double> rho2;
  lucia_cuda::DeviceBuffer<double> rho2s;
  lucia_cuda::DeviceBuffer<double> rho2a;
  RouteSlot route[2];
  unsigned int next_route = 0;
  lucia_cuda::PinnedBuffer<double> staging;
  bool device_ready = false;
  bool session_active = false;
  bool session_resident = false;
  double *session_x_host = nullptr;
  const double *session_sb_host = nullptr;
  const double *session_cb_host = nullptr;
  std::size_t session_x_count = 0;
  std::size_t session_sb_count = 0;
  std::size_t session_cb_count = 0;
  bool density_active = false;
  bool density_dirty = false;
  bool density_ipack = false;
  std::size_t density_nacob = 0;
  std::size_t density_count = 0;
  double *density_rho2_host = nullptr;
  double *density_rho2s_host = nullptr;
  double *density_rho2a_host = nullptr;
  std::size_t max_grid_x = 0;
  std::size_t max_grid_y = 0;
  std::size_t max_grid_z = 0;

  bool session_matches(double *x_host, const double *sb_host, const double *cb_host, std::size_t x_count, std::size_t sb_count,
                       std::size_t cb_count) const noexcept {
    return session_active && session_x_host == x_host && session_sb_host == sb_host && session_cb_host == cb_host
           && session_x_count == x_count && session_sb_count == sb_count && session_cb_count == cb_count;
  }

  bool begin_session(double *x_host, const double *sb_host, const double *cb_host, std::size_t x_count, std::size_t sb_count,
                     std::size_t cb_count) noexcept {
    if (session_active) {
      return false;
    }
    session_active = true;
    session_resident = false;
    session_x_host = x_host;
    session_sb_host = sb_host;
    session_cb_host = cb_host;
    session_x_count = x_count;
    session_sb_count = sb_count;
    session_cb_count = cb_count;
    return true;
  }

  bool flush_session() noexcept {
    if (!session_active || !session_resident) {
      return true;
    }
    invalidate_residency();
    std::size_t x_bytes = 0;
    if (!lucia_cuda::checked_mul(session_x_count, sizeof(double), &x_bytes) || !staging.reserve(session_x_count)
        || cudaMemcpy(staging.data, x.data, x_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) {
      return false;
    }
    std::memcpy(session_x_host, staging.data, x_bytes);
    return true;
  }

  void invalidate_residency() noexcept {
    session_resident = false;
  }

  std::int64_t fallback_status() noexcept {
    if (!flush_session() || !disable_density()) {
      return -1;
    }
    return 0;
  }

  int begin_density(double *rho2_host, double *rho2s_host, double *rho2a_host, std::size_t nacob, bool ipack,
                    std::size_t density_count) noexcept {
    if (density_active) {
      return 0;
    }
    density_ipack = ipack;
    density_nacob = nacob;
    this->density_count = density_count;
    density_rho2_host = rho2_host;
    density_rho2s_host = rho2s_host;
    density_rho2a_host = rho2a_host;
    const bool copied = ipack ? rho2s.copy_from(rho2s_host, density_count) && rho2a.copy_from(rho2a_host, density_count)
                              : rho2.copy_from(rho2_host, density_count);
    if (!copied) {
      clear_density();
      return 0;
    }
    density_active = true;
    density_dirty = false;
    return 1;
  }

  bool flush_density() noexcept {
    if (!density_active || !density_dirty) {
      return true;
    }
    std::size_t bytes = 0;
    if (!lucia_cuda::checked_mul(density_count, sizeof(double), &bytes)) {
      return false;
    }
    if (density_ipack) {
      if (cudaMemcpy(density_rho2s_host, rho2s.data, bytes, cudaMemcpyDeviceToHost) != cudaSuccess
          || cudaMemcpy(density_rho2a_host, rho2a.data, bytes, cudaMemcpyDeviceToHost) != cudaSuccess) {
        return false;
      }
    } else if (cudaMemcpy(density_rho2_host, rho2.data, bytes, cudaMemcpyDeviceToHost) != cudaSuccess) {
      return false;
    }
    density_dirty = false;
    return true;
  }

  void clear_density() noexcept {
    density_active = false;
    density_dirty = false;
    density_ipack = false;
    density_nacob = 0;
    density_count = 0;
    density_rho2_host = nullptr;
    density_rho2s_host = nullptr;
    density_rho2a_host = nullptr;
  }

  bool disable_density() noexcept {
    const bool success = flush_density();
    if (success) {
      clear_density();
    }
    return success;
  }

  bool end_density() noexcept {
    return disable_density();
  }

  void clear_session() noexcept {
    session_active = false;
    session_resident = false;
    session_x_host = nullptr;
    session_sb_host = nullptr;
    session_cb_host = nullptr;
    session_x_count = 0;
    session_sb_count = 0;
    session_cb_count = 0;
  }

  bool end_session() noexcept {
    const bool success = flush_session();
    clear_session();
    return success;
  }

  bool supports_grid(std::size_t x_count, std::size_t y_count, std::size_t z_count) noexcept {
    if (!device_ready) {
      int device = 0;
      cudaDeviceProp properties{};
      if (cudaGetDevice(&device) != cudaSuccess || cudaGetDeviceProperties(&properties, device) != cudaSuccess
          || properties.maxGridSize[0] <= 0 || properties.maxGridSize[1] <= 0 || properties.maxGridSize[2] <= 0) {
        return false;
      }
      max_grid_x = static_cast<std::size_t>(properties.maxGridSize[0]);
      max_grid_y = static_cast<std::size_t>(properties.maxGridSize[1]);
      max_grid_z = static_cast<std::size_t>(properties.maxGridSize[2]);
      device_ready = true;
    }
    const std::size_t max_unsigned = static_cast<std::size_t>((std::numeric_limits<unsigned int>::max)());
    return x_count > 0 && y_count > 0 && z_count > 0 && x_count <= max_grid_x && y_count <= max_grid_y && z_count <= max_grid_z
           && x_count <= max_unsigned && y_count <= max_unsigned && z_count <= max_unsigned;
  }

  bool release() noexcept {
    bool success = end_session();
    success = end_density() && success;
    if (!x.release())
      success = false;
    if (!sb.release())
      success = false;
    if (!cb.release())
      success = false;
    if (!rho2.release())
      success = false;
    if (!rho2s.release())
      success = false;
    if (!rho2a.release())
      success = false;
    if (!route[0].release())
      success = false;
    if (!route[1].release())
      success = false;
    next_route = 0;
    if (!staging.release())
      success = false;
    device_ready = false;
    max_grid_x = 0;
    max_grid_y = 0;
    max_grid_z = 0;
    return success;
  }
};

__device__ double atomic_add(double *address, double value) {
#if __CUDA_ARCH__ >= 600
  return atomicAdd(address, value);
#else
  auto bits = reinterpret_cast<unsigned long long *>(address);
  unsigned long long old = *bits;
  unsigned long long assumed = 0;
  do {
    assumed = old;
    old = atomicCAS(bits, assumed, __double_as_longlong(value + __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
#endif
}

__device__ std::size_t triangular_count(std::size_t n) {
  return (n & 1U) == 0U ? (n / 2) * (n + 1) : n * ((n + 1) / 2);
}

__device__ void decode_triangular(std::size_t index, std::size_t extent, std::size_t *first, std::size_t *second) {
  std::size_t lower = 1;
  std::size_t upper = extent;
  while (lower < upper) {
    const std::size_t middle = lower + (upper - lower) / 2;
    if (triangular_count(middle) > index) {
      upper = middle;
    } else {
      lower = middle + 1;
    }
  }
  *first = lower;
  *second = index - triangular_count(lower - 1) + 1;
}

__device__ std::size_t triangular_index(std::size_t first, std::size_t second) {
  const std::size_t high = first > second ? first : second;
  const std::size_t low = first > second ? second : first;
  return high * (high - 1) / 2 + low;
}

__device__ std::size_t packed_pair_index(std::size_t first, std::size_t second) {
  return first * (first - 1) / 2 + second;
}

__device__ void scatter_density(double *rho2, double *rho2s, double *rho2a, std::size_t nacob, bool ipack, std::size_t i,
                                std::size_t j, std::size_t k, std::size_t l, double value, double sign) {
  const std::size_t ij = (j - 1) * nacob + i;
  const std::size_t kl = (l - 1) * nacob + k;
  if (!ipack) {
    if (ij >= kl) {
      atomic_add(rho2 + triangular_index(ij, kl) - 1, -sign * value);
    }
    return;
  }

  if (k < l) {
    return;
  }
  const std::size_t ij_pack = packed_pair_index(i, j);
  const std::size_t ji_pack = packed_pair_index(j, i);
  const std::size_t kl_pack = packed_pair_index(k, l);
  const double term = (k == l ? 0.25 : 0.5) * sign * value;
  if (i >= j && ij_pack >= kl_pack) {
    const std::size_t target = packed_pair_index(ij_pack, kl_pack) - 1;
    atomic_add(rho2s + target, -term);
    atomic_add(rho2a + target, -term);
  }
  if (j >= i && ji_pack >= kl_pack) {
    const std::size_t target = packed_pair_index(ji_pack, kl_pack) - 1;
    atomic_add(rho2s + target, -term);
    atomic_add(rho2a + target, term);
  }
}

__global__ void scale_cuda_kernel(double *x, std::size_t count, double factor) {
  const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (index < count) {
    x[index] *= factor;
  }
}

__global__ void gsbbd2a_density_kernel(const double *x, double *rho2, double *rho2s, double *rho2a, std::size_t nik,
                                       std::size_t njl, std::size_t ni, std::size_t nk, std::size_t nj, std::size_t nl,
                                       std::size_t ioff, std::size_t joff, std::size_t koff, std::size_t loff, std::size_t nacob,
                                       bool ik_tri, bool jl_tri, bool ipack) {
  const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  const std::size_t count = nik * njl;
  if (index >= count) {
    return;
  }

  std::size_t i = 0;
  std::size_t k = 0;
  std::size_t j = 0;
  std::size_t l = 0;
  const std::size_t ik_index = index % nik;
  const std::size_t jl_index = index / nik;
  if (ik_tri) {
    decode_triangular(ik_index, ni, &i, &k);
  } else {
    i = ik_index % ni + 1;
    k = ik_index / ni + 1;
  }
  if (jl_tri) {
    decode_triangular(jl_index, nj, &j, &l);
  } else {
    j = jl_index % nj + 1;
    l = jl_index / nj + 1;
  }

  const std::size_t i_variants = ik_tri && i != k ? 2 : 1;
  const std::size_t j_variants = jl_tri && j != l ? 2 : 1;
  const std::size_t i_global = i + ioff - 1;
  const std::size_t k_global = k + koff - 1;
  const std::size_t j_global = j + joff - 1;
  const std::size_t l_global = l + loff - 1;
  for (std::size_t i_variant = 0; i_variant < i_variants; ++i_variant) {
    const std::size_t output_i = i_variant == 0 ? i_global : k_global;
    const std::size_t output_k = i_variant == 0 ? k_global : i_global;
    const double sign_ik = i_variant == 0 ? 1.0 : -1.0;
    for (std::size_t j_variant = 0; j_variant < j_variants; ++j_variant) {
      const std::size_t output_j = j_variant == 0 ? j_global : l_global;
      const std::size_t output_l = j_variant == 0 ? l_global : j_global;
      const double sign = sign_ik * (j_variant == 0 ? 1.0 : -1.0);
      scatter_density(rho2, rho2s, rho2a, nacob, ipack, output_i, output_j, output_k, output_l, x[index], sign);
      if (!ik_tri) {
        scatter_density(rho2, rho2s, rho2a, nacob, ipack, output_k, output_j, output_i, output_l, x[index], -sign);
      }
      if (!jl_tri) {
        scatter_density(rho2, rho2s, rho2a, nacob, ipack, output_i, output_l, output_k, output_j, x[index], -sign);
      }
      if (!ik_tri && !jl_tri) {
        scatter_density(rho2, rho2s, rho2a, nacob, ipack, output_k, output_l, output_i, output_j, x[index], sign);
      }
    }
  }
}

__global__ void gsbbd2a_cuda_kernel(double *x, const double *sb, const double *cb, const std::int64_t *i1, const double *xi1s,
                                    const std::int64_t *i2, const double *xi2s, const std::int64_t *ik_map_col,
                                    const std::int64_t *jl_map_col, std::size_t nrow, std::size_t ibot0, std::size_t nibtc,
                                    std::size_t maxk, std::size_t nik, std::size_t njl, std::size_t reduction_count,
                                    std::size_t reduction_per_split) {
  __shared__ double a_tile[tile_size][tile_size + 1];
  __shared__ double b_tile[tile_size][tile_size + 1];

  const std::size_t ik = static_cast<std::size_t>(blockIdx.x) * tile_size + threadIdx.y;
  const std::size_t jl = static_cast<std::size_t>(blockIdx.y) * tile_size + threadIdx.x;
  double sum = 0.0;

  std::size_t base = static_cast<std::size_t>(blockIdx.z) * reduction_per_split;
  if (base >= reduction_count) {
    return;
  }
  const std::size_t split_size = reduction_per_split < reduction_count - base ? reduction_per_split : reduction_count - base;
  const std::size_t reduction_end = base + split_size;
  while (base < reduction_end) {
    const std::size_t remaining = reduction_end - base;
    const bool a_active = static_cast<std::size_t>(threadIdx.x) < remaining;
    double a_value = 0.0;
    if (a_active && ik < nik && ik_map_col[ik] >= 0) {
      const std::size_t a_k = base + static_cast<std::size_t>(threadIdx.x);
      const std::size_t q = a_k / nibtc;
      const std::size_t r = a_k % nibtc;
      const std::int64_t map_col = ik_map_col[ik];
      const std::size_t map_offset = q + static_cast<std::size_t>(map_col) * maxk;
      const std::int64_t source = i1[map_offset];
      if (source != 0) {
        const std::size_t source_offset = ibot0 + r + static_cast<std::size_t>(source - 1) * nrow;
        a_value = xi1s[map_offset] * sb[source_offset];
      }
    }
    a_tile[threadIdx.y][threadIdx.x] = a_value;

    const bool b_active = static_cast<std::size_t>(threadIdx.y) < remaining;
    double b_value = 0.0;
    if (b_active && jl < njl && jl_map_col[jl] >= 0) {
      const std::size_t b_k = base + static_cast<std::size_t>(threadIdx.y);
      const std::size_t q = b_k / nibtc;
      const std::size_t r = b_k % nibtc;
      const std::int64_t map_col = jl_map_col[jl];
      const std::size_t map_offset = q + static_cast<std::size_t>(map_col) * maxk;
      const std::int64_t source = i2[map_offset];
      if (source != 0) {
        const std::size_t source_offset = ibot0 + r + static_cast<std::size_t>(source - 1) * nrow;
        b_value = xi2s[map_offset] * cb[source_offset];
      }
    }
    b_tile[threadIdx.y][threadIdx.x] = b_value;

    __syncthreads();
    if (ik < nik && jl < njl) {
      for (std::size_t k = 0; k < tile_size; ++k) {
        sum += a_tile[threadIdx.y][k] * b_tile[k][threadIdx.x];
      }
    }
    __syncthreads();

    base += tile_size;
  }

  if (ik < nik && jl < njl) {
    const std::size_t output = ik + jl * nik;
    atomic_add(x + output, -sum);
  }
}

Workspace workspace;

std::int64_t block_fallback() noexcept {
  const std::int64_t status = workspace.fallback_status();
  if (status != 0) {
    return -1;
  }
  workspace.clear_session();
  return 0;
}

std::int64_t finish_block(std::int64_t ni, std::int64_t ioff, std::int64_t nj, std::int64_t joff, std::int64_t nk,
                          std::int64_t koff, std::int64_t nl, std::int64_t loff, std::int64_t nacob, std::int64_t ipack) noexcept {
  if (!workspace.session_active) {
    return workspace.disable_density() ? 0 : -1;
  }
  if (!workspace.session_resident) {
    if (!workspace.disable_density()) {
      return -1;
    }
    workspace.clear_session();
    return 0;
  }
  if (!workspace.density_active) {
    return workspace.end_session() ? 0 : -1;
  }

  std::size_t ni_size = 0, ioff_size = 0, nj_size = 0, joff_size = 0;
  std::size_t nk_size = 0, koff_size = 0, nl_size = 0, loff_size = 0;
  std::size_t nacob_size = 0;
  if (!lucia_cuda::positive_size(ni, &ni_size) || !lucia_cuda::positive_size(ioff, &ioff_size)
      || !lucia_cuda::positive_size(nj, &nj_size) || !lucia_cuda::positive_size(joff, &joff_size)
      || !lucia_cuda::positive_size(nk, &nk_size) || !lucia_cuda::positive_size(koff, &koff_size)
      || !lucia_cuda::positive_size(nl, &nl_size) || !lucia_cuda::positive_size(loff, &loff_size)
      || !lucia_cuda::positive_size(nacob, &nacob_size) || (ipack != 0 && ipack != 1) || nacob_size != workspace.density_nacob
      || (ipack == 1) != workspace.density_ipack || !valid_extent(ioff_size, ni_size, nacob_size)
      || !valid_extent(koff_size, nk_size, nacob_size) || !valid_extent(joff_size, nj_size, nacob_size)
      || !valid_extent(loff_size, nl_size, nacob_size)) {
    return block_fallback();
  }

  const bool ik_tri = ioff_size == koff_size;
  const bool jl_tri = joff_size == loff_size;
  if ((ik_tri && ni_size != nk_size) || (jl_tri && nj_size != nl_size)) {
    return block_fallback();
  }
  std::size_t nik = 0, njl = 0, x_count = 0;
  if ((ik_tri ? !checked_triangular(ni_size, &nik) : !lucia_cuda::checked_mul(ni_size, nk_size, &nik))
      || (jl_tri ? !checked_triangular(nj_size, &njl) : !lucia_cuda::checked_mul(nj_size, nl_size, &njl))
      || !lucia_cuda::checked_mul(nik, njl, &x_count)) {
    return block_fallback();
  }
  std::size_t grid_x = 0;
  if (!checked_ceil_div(x_count, scale_block_size, &grid_x) || !workspace.supports_grid(grid_x, 1, 1)) {
    return block_fallback();
  }

  gsbbd2a_density_kernel<<<static_cast<unsigned int>(grid_x), static_cast<unsigned int>(scale_block_size)>>>(
      workspace.x.data, workspace.rho2.data, workspace.rho2s.data, workspace.rho2a.data, nik, njl, ni_size, nk_size, nj_size,
      nl_size, ioff_size, joff_size, koff_size, loff_size, nacob_size, ik_tri, jl_tri, ipack == 1);
  if (cudaGetLastError() != cudaSuccess) {
    workspace.invalidate_residency();
    return -1;
  }
  workspace.clear_session();
  workspace.density_dirty = true;
  return 1;
}

} // namespace

extern "C" int64_t lucia_gsbbd2a_cuda_begin(double *x, const double *sb, const double *cb, int64_t nrow, int64_t nsb, int64_t ncb,
                                            int64_t xcount) {
  if (x == nullptr || sb == nullptr || cb == nullptr || workspace.session_active) {
    if (!workspace.disable_density()) {
      return -1;
    }
    return 0;
  }
  std::size_t nrow_size = 0;
  std::size_t nsb_size = 0;
  std::size_t ncb_size = 0;
  std::size_t x_count = 0;
  std::size_t sb_count = 0;
  std::size_t cb_count = 0;
  if (!lucia_cuda::positive_size(nrow, &nrow_size) || !lucia_cuda::positive_size(nsb, &nsb_size)
      || !lucia_cuda::positive_size(ncb, &ncb_size) || !lucia_cuda::positive_size(xcount, &x_count)
      || !lucia_cuda::checked_mul(nrow_size, nsb_size, &sb_count) || !lucia_cuda::checked_mul(nrow_size, ncb_size, &cb_count)
      || !workspace.begin_session(x, sb, cb, x_count, sb_count, cb_count)) {
    if (!workspace.disable_density()) {
      return -1;
    }
    return 0;
  }
  return 1;
}

extern "C" int64_t lucia_gsbbd2a_cuda_density_begin(double *rho2, double *rho2s, double *rho2a, int64_t nacob, int64_t ipack) {
  if (rho2 == nullptr || (ipack != 0 && ipack != 1) || (ipack == 1 && (rho2s == nullptr || rho2a == nullptr))) {
    return 0;
  }
  std::size_t nacob_size = 0;
  std::size_t density_count = 0;
  std::size_t packed_count = 0;
  if (!lucia_cuda::positive_size(nacob, &nacob_size)
      || (ipack == 0
              ? !lucia_cuda::checked_mul(nacob_size, nacob_size, &packed_count) || !checked_triangular(packed_count, &density_count)
              : !checked_triangular(nacob_size, &packed_count) || !checked_triangular(packed_count, &density_count))) {
    return 0;
  }
  return workspace.begin_density(rho2, rho2s, rho2a, nacob_size, ipack == 1, density_count);
}

extern "C" int64_t lucia_gsbbd2a_cuda_density_end() {
  return workspace.end_density() ? 1 : -1;
}

extern "C" int64_t lucia_gsbbd2a_cuda_block_end(int64_t ni, int64_t ioff, int64_t nj, int64_t joff, int64_t nk, int64_t koff,
                                                int64_t nl, int64_t loff, int64_t nacob, int64_t ipack) {
  return finish_block(ni, ioff, nj, joff, nk, koff, nl, loff, nacob, ipack);
}

extern "C" int64_t lucia_gsbbd2a_cuda_end() {
  return workspace.end_session() ? 1 : -1;
}

extern "C" int64_t lucia_gsbbd2a_cuda_route(double *x, const double *sb, const double *cb, const int64_t *i1, const double *xi1s,
                                            const int64_t *i2, const double *xi2s, int64_t nrow, int64_t nsb, int64_t ncb,
                                            int64_t ibot, int64_t nibtc, int64_t maxk, int64_t nkbtc, int64_t ni, int64_t nk,
                                            int64_t nj, int64_t nl, int64_t nik, int64_t njl, int64_t iksm, int64_t jlsm,
                                            double factor) {
  if (x == nullptr || sb == nullptr || cb == nullptr || i1 == nullptr || xi1s == nullptr || i2 == nullptr || xi2s == nullptr
      || !std::isfinite(factor)) {
    return workspace.fallback_status();
  }

  std::size_t nrow_size = 0;
  std::size_t nsb_size = 0;
  std::size_t ncb_size = 0;
  std::size_t ibot_size = 0;
  std::size_t nibtc_size = 0;
  std::size_t maxk_size = 0;
  std::size_t nkbtc_size = 0;
  std::size_t ni_size = 0;
  std::size_t nk_size = 0;
  std::size_t nj_size = 0;
  std::size_t nl_size = 0;
  std::size_t nik_size = 0;
  std::size_t njl_size = 0;
  if (!lucia_cuda::positive_size(nrow, &nrow_size) || !lucia_cuda::positive_size(nsb, &nsb_size)
      || !lucia_cuda::positive_size(ncb, &ncb_size) || !lucia_cuda::positive_size(ibot, &ibot_size)
      || !lucia_cuda::positive_size(nibtc, &nibtc_size) || !lucia_cuda::positive_size(maxk, &maxk_size)
      || !lucia_cuda::positive_size(nkbtc, &nkbtc_size) || !lucia_cuda::positive_size(ni, &ni_size)
      || !lucia_cuda::positive_size(nk, &nk_size) || !lucia_cuda::positive_size(nj, &nj_size)
      || !lucia_cuda::positive_size(nl, &nl_size) || !lucia_cuda::positive_size(nik, &nik_size)
      || !lucia_cuda::positive_size(njl, &njl_size)) {
    return workspace.fallback_status();
  }
  if (nkbtc_size > maxk_size || iksm < 0 || iksm > 1 || jlsm < 0 || jlsm > 1) {
    return workspace.fallback_status();
  }

  const std::size_t ibot0 = ibot_size - 1;
  if (ibot0 >= nrow_size || nibtc_size > nrow_size - ibot0) {
    return workspace.fallback_status();
  }

  std::size_t full_ik = 0;
  std::size_t full_jl = 0;
  std::size_t x_block_count = 0;
  std::size_t x_count = 0;
  std::size_t sb_count = 0;
  std::size_t cb_count = 0;
  std::size_t i1_count = 0;
  std::size_t xi1s_count = 0;
  std::size_t i2_count = 0;
  std::size_t xi2s_count = 0;
  std::size_t reduction_count = 0;
  if (!lucia_cuda::checked_mul(ni_size, nk_size, &full_ik) || !lucia_cuda::checked_mul(nj_size, nl_size, &full_jl)
      || !lucia_cuda::checked_mul(ni_size, nj_size, &x_block_count)
      || !lucia_cuda::checked_mul(x_block_count, nk_size, &x_block_count)
      || !lucia_cuda::checked_mul(x_block_count, nl_size, &x_block_count) || !lucia_cuda::checked_mul(nik_size, njl_size, &x_count)
      || !lucia_cuda::checked_mul(nrow_size, nsb_size, &sb_count) || !lucia_cuda::checked_mul(nrow_size, ncb_size, &cb_count)
      || !lucia_cuda::checked_mul(maxk_size, full_ik, &i1_count) || !lucia_cuda::checked_mul(maxk_size, full_ik, &xi1s_count)
      || !lucia_cuda::checked_mul(maxk_size, full_jl, &i2_count) || !lucia_cuda::checked_mul(maxk_size, full_jl, &xi2s_count)
      || !lucia_cuda::checked_mul(nibtc_size, nkbtc_size, &reduction_count)) {
    return workspace.fallback_status();
  }
  std::size_t x_bytes = 0;
  if (!lucia_cuda::checked_mul(x_count, sizeof(double), &x_bytes)) {
    return workspace.fallback_status();
  }

  const bool session_match = workspace.session_matches(x, sb, cb, x_block_count, sb_count, cb_count);
  if (workspace.session_active && !session_match && workspace.fallback_status() != 0) {
    return -1;
  }
  const bool shared_match = lucia_sigma_cuda_blocks::sigma_blocks_match(const_cast<double *>(sb), cb, sb_count, cb_count);

  std::size_t grid_x = 0;
  std::size_t grid_y = 0;
  std::size_t scale_grid_x = 0;
  std::size_t split_count = 0;
  std::size_t reduction_per_split = 0;
  if (!checked_ceil_div(nik_size, tile_size, &grid_x) || !checked_ceil_div(njl_size, tile_size, &grid_y)
      || !checked_ceil_div(x_count, scale_block_size, &scale_grid_x)
      || !checked_ceil_div(reduction_count, reduction_target_per_split, &split_count)) {
    return workspace.fallback_status();
  }
  if (split_count > max_split_count) {
    split_count = max_split_count;
  }
  if (!checked_ceil_div(reduction_count, split_count, &reduction_per_split)) {
    return workspace.fallback_status();
  }

  if (!validate_map_columns(i1, iksm, ni_size, nk_size, nik_size, maxk_size, nkbtc_size, nsb)
      || !validate_map_columns(i2, jlsm, nj_size, nl_size, njl_size, maxk_size, nkbtc_size, ncb)) {
    return workspace.fallback_status();
  }

  const bool upload_blocks = session_match && !workspace.session_resident;
  if (!workspace.supports_grid(grid_x, grid_y, split_count) || !workspace.supports_grid(scale_grid_x, 1, 1)) {
    return workspace.fallback_status();
  }

  std::size_t i1_bytes = 0, xi1s_bytes = 0, i2_bytes = 0, xi2s_bytes = 0;
  std::size_t ik_map_col_bytes = 0, jl_map_col_bytes = 0;
  if (!lucia_cuda::checked_mul(i1_count, sizeof(std::int64_t), &i1_bytes)
      || !lucia_cuda::checked_mul(xi1s_count, sizeof(double), &xi1s_bytes)
      || !lucia_cuda::checked_mul(i2_count, sizeof(std::int64_t), &i2_bytes)
      || !lucia_cuda::checked_mul(xi2s_count, sizeof(double), &xi2s_bytes)
      || !lucia_cuda::checked_mul(nik_size, sizeof(std::int64_t), &ik_map_col_bytes)
      || !lucia_cuda::checked_mul(njl_size, sizeof(std::int64_t), &jl_map_col_bytes)) {
    return workspace.fallback_status();
  }
  std::size_t i1_offset = 0, xi1s_offset = 0, i2_offset = 0, xi2s_offset = 0;
  std::size_t ik_map_col_offset = 0, jl_map_col_offset = 0;
  std::size_t payload_bytes = 0;
  if (!append_section(i1_bytes, &i1_offset, &payload_bytes) || !append_section(xi1s_bytes, &xi1s_offset, &payload_bytes)
      || !append_section(i2_bytes, &i2_offset, &payload_bytes) || !append_section(xi2s_bytes, &xi2s_offset, &payload_bytes)
      || !append_section(ik_map_col_bytes, &ik_map_col_offset, &payload_bytes)
      || !append_section(jl_map_col_bytes, &jl_map_col_offset, &payload_bytes)) {
    return workspace.fallback_status();
  }

  RouteSlot &slot = workspace.route[workspace.next_route];
  if (!slot.prepare() || !slot.host.reserve(payload_bytes)) {
    return workspace.fallback_status();
  }
  unsigned char *payload_host = slot.host.data;
  std::memcpy(payload_host + i1_offset, i1, i1_bytes);
  std::memcpy(payload_host + xi1s_offset, xi1s, xi1s_bytes);
  std::memcpy(payload_host + i2_offset, i2, i2_bytes);
  std::memcpy(payload_host + xi2s_offset, xi2s, xi2s_bytes);
  std::int64_t *ik_map_col_host = reinterpret_cast<std::int64_t *>(payload_host + ik_map_col_offset);
  std::int64_t *jl_map_col_host = reinterpret_cast<std::int64_t *>(payload_host + jl_map_col_offset);
  if (!fill_pair_columns(iksm, ni_size, nk_size, nik_size, ik_map_col_host)
      || !fill_pair_columns(jlsm, nj_size, nl_size, njl_size, jl_map_col_host)) {
    return workspace.fallback_status();
  }

  double *device_sb = workspace.sb.data;
  const double *device_cb = workspace.cb.data;
  if (shared_match
      && !lucia_sigma_cuda_blocks::sigma_blocks_acquire(const_cast<double *>(sb), cb, sb_count, cb_count, &device_sb, &device_cb)) {
    return workspace.fallback_status();
  }
  if ((!session_match || upload_blocks)
      && (!workspace.staging.reserve(session_match ? x_block_count : x_count)
          || !workspace.x.copy_from(x, session_match ? x_block_count : x_count))) {
    return workspace.fallback_status();
  }
  if (!shared_match && (!session_match || upload_blocks)
      && (!workspace.sb.copy_from(sb, sb_count) || !workspace.cb.copy_from(cb, cb_count))) {
    return workspace.fallback_status();
  }
  if (!shared_match) {
    device_sb = workspace.sb.data;
    device_cb = workspace.cb.data;
  }
  if (!slot.upload(payload_bytes)) {
    return workspace.fallback_status();
  }
  workspace.next_route = (workspace.next_route + 1U) % 2U;

  const unsigned char *payload_device = slot.device.data;
  const std::int64_t *i1_device = reinterpret_cast<const std::int64_t *>(payload_device + i1_offset);
  const double *xi1s_device = reinterpret_cast<const double *>(payload_device + xi1s_offset);
  const std::int64_t *i2_device = reinterpret_cast<const std::int64_t *>(payload_device + i2_offset);
  const double *xi2s_device = reinterpret_cast<const double *>(payload_device + xi2s_offset);
  const std::int64_t *ik_map_col_device = reinterpret_cast<const std::int64_t *>(payload_device + ik_map_col_offset);
  const std::int64_t *jl_map_col_device = reinterpret_cast<const std::int64_t *>(payload_device + jl_map_col_offset);

  const dim3 block(static_cast<unsigned int>(tile_size), static_cast<unsigned int>(tile_size), 1);
  const dim3 grid(static_cast<unsigned int>(grid_x), static_cast<unsigned int>(grid_y), static_cast<unsigned int>(split_count));
  scale_cuda_kernel<<<static_cast<unsigned int>(scale_grid_x), static_cast<unsigned int>(scale_block_size)>>>(workspace.x.data,
                                                                                                              x_count, factor);
  gsbbd2a_cuda_kernel<<<grid, block>>>(workspace.x.data, device_sb, device_cb, i1_device, xi1s_device, i2_device, xi2s_device,
                                       ik_map_col_device, jl_map_col_device, nrow_size, ibot0, nibtc_size, maxk_size, nik_size,
                                       njl_size, reduction_count, reduction_per_split);

  bool success = cudaGetLastError() == cudaSuccess;
  if (!session_match && cudaDeviceSynchronize() != cudaSuccess) {
    success = false;
  }
  if (!success) {
    workspace.invalidate_residency();
    return -1;
  }
  if (session_match) {
    workspace.session_resident = true;
    return 1;
  }

  if (cudaMemcpy(workspace.staging.data, workspace.x.data, x_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) {
    return -1;
  }

  std::memcpy(x, workspace.staging.data, x_bytes);
  return 1;
}

extern "C" void lucia_gsbbd2a_cuda_release() {
  workspace.release();
}
