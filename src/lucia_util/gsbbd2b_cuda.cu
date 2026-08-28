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
constexpr std::size_t maps_per_split = 8;
constexpr std::size_t minimum_active_work = 10000000;

bool checked_add(std::size_t left, std::size_t right, std::size_t *result) noexcept
{
  if (left > (std::numeric_limits<std::size_t>::max)() - right) {
    return false;
  }
  *result = left + right;
  return true;
}
bool checked_ntri(std::size_t n, std::size_t *result) noexcept
{
  std::size_t next = 0;
  if (!checked_add(n, 1, &next)) {
    return false;
  }
  if ((n & 1U) == 0U) {
    return lucia_cuda::checked_mul(n / 2, next, result);
  }
  return lucia_cuda::checked_mul(n, next / 2, result);
}

bool valid_extent(std::size_t offset, std::size_t extent, std::size_t limit) noexcept
{
  return offset > 0 && extent > 0 && offset <= limit && extent <= limit - (offset - 1);
}

template <typename T>
int upload_device(lucia_cuda::DeviceBuffer<T> &device, const T *host, std::size_t count) noexcept
{
  if (count == 0) {
    return 1;
  }
  if (host == nullptr || !device.reserve(count)) {
    return 0;
  }
  std::size_t bytes = 0;
  if (!lucia_cuda::checked_mul(count, sizeof(T), &bytes)) {
    return 0;
  }
  return cudaMemcpy(device.data, host, bytes, cudaMemcpyHostToDevice) == cudaSuccess ? 1 : -1;
}

struct Workspace {
  lucia_cuda::DeviceBuffer<double> x;
  lucia_cuda::DeviceBuffer<double> sb;
  lucia_cuda::DeviceBuffer<double> cb;
  lucia_cuda::DeviceBuffer<double> xi1;
  lucia_cuda::DeviceBuffer<double> xi3;
  lucia_cuda::DeviceBuffer<double> xi4;
  lucia_cuda::DeviceBuffer<double> xi2;
  lucia_cuda::DeviceBuffer<std::int64_t> i1;
  lucia_cuda::DeviceBuffer<std::int64_t> i3;
  lucia_cuda::DeviceBuffer<std::int64_t> i4;
  lucia_cuda::DeviceBuffer<std::int64_t> i2;
  lucia_cuda::DeviceBuffer<double> rho2;
  lucia_cuda::DeviceBuffer<double> rho2s;
  lucia_cuda::DeviceBuffer<double> rho2a;
  lucia_cuda::DeviceBuffer<double> s2_delta;
  lucia_cuda::PinnedBuffer<double> staging;
  bool device_ready = false;
  bool session_active = false;
  bool session_resident = false;
  const double *session_sb_host = nullptr;
  const double *session_cb_host = nullptr;
  std::size_t session_sb_count = 0;
  std::size_t session_cb_count = 0;
  bool density_active = false;
  bool density_uploaded = false;
  bool density_dirty = false;
  bool density_ipack = false;
  std::size_t density_norb = 0;
  std::size_t density_count = 0;
  double *density_rho2_host = nullptr;
  double *density_rho2s_host = nullptr;
  double *density_rho2a_host = nullptr;
  double *density_s2_host = nullptr;
  bool map_session_active = false;
  bool map_session_resident = false;
  const std::int64_t *map_i1_host = nullptr;
  const double *map_xi1_host = nullptr;
  const std::int64_t *map_i3_host = nullptr;
  const double *map_xi3_host = nullptr;
  std::size_t map_i1_count = 0;
  std::size_t map_i3_count = 0;
  std::size_t max_grid_x = 0;
  std::size_t max_grid_y = 0;
  std::size_t max_grid_z = 0;

  bool session_matches(const double *sb_host, const double *cb_host,
                       std::size_t sb_count, std::size_t cb_count) const noexcept
  {
    return session_active && session_sb_host == sb_host && session_cb_host == cb_host &&
           session_sb_count == sb_count && session_cb_count == cb_count;
  }

  bool begin_session(const double *sb_host, const double *cb_host,
                     std::size_t sb_count, std::size_t cb_count) noexcept
  {
    if (session_active) {
      return false;
    }
    session_active = true;
    session_resident = false;
    session_sb_host = sb_host;
    session_cb_host = cb_host;
    session_sb_count = sb_count;
    session_cb_count = cb_count;
    return true;
  }

  void clear_session() noexcept
  {
    session_active = false;
    session_resident = false;
    session_sb_host = nullptr;
    session_cb_host = nullptr;
    session_sb_count = 0;
    session_cb_count = 0;
  }

  bool end_session() noexcept
  {
    const bool success = flush_density();
    clear_density_state();
    clear_session();
    return success;
  }

  bool begin_density_session(double *rho2_host, double *rho2s_host, double *rho2a_host,
                             double *s2_host, std::size_t norb, bool ipack,
                             std::size_t density_count) noexcept
  {
    density_active = true;
    density_uploaded = false;
    density_dirty = false;
    density_ipack = ipack;
    density_norb = norb;
    this->density_count = density_count;
    density_rho2_host = rho2_host;
    density_rho2s_host = rho2s_host;
    density_rho2a_host = rho2a_host;
    density_s2_host = s2_host;
    const bool copied = ipack
                            ? rho2s.copy_from(rho2s_host, density_count) &&
                                  rho2a.copy_from(rho2a_host, density_count)
                            : rho2.copy_from(rho2_host, density_count);
    if (!copied || !s2_delta.reserve(1) ||
        cudaMemset(s2_delta.data, 0, sizeof(double)) != cudaSuccess) {
      clear_density_state();
      return false;
    }
    density_uploaded = true;
    return true;
  }

  bool flush_density() noexcept
  {
    if (!density_active || !density_dirty) {
      return true;
    }
    std::size_t density_bytes = 0;
    if (!lucia_cuda::checked_mul(density_count, sizeof(double), &density_bytes) ||
        density_s2_host == nullptr) {
      return false;
    }
    if (density_ipack) {
      if (cudaMemcpy(density_rho2s_host, rho2s.data, density_bytes,
                     cudaMemcpyDeviceToHost) != cudaSuccess ||
          cudaMemcpy(density_rho2a_host, rho2a.data, density_bytes,
                     cudaMemcpyDeviceToHost) != cudaSuccess) {
        return false;
      }
    } else if (cudaMemcpy(density_rho2_host, rho2.data, density_bytes,
                          cudaMemcpyDeviceToHost) != cudaSuccess) {
      return false;
    }
    double delta = 0.0;
    if (cudaMemcpy(&delta, s2_delta.data, sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess) {
      return false;
    }
    *density_s2_host += delta;
    density_dirty = false;
    return true;
  }

  void clear_density_state() noexcept
  {
    density_active = false;
    density_uploaded = false;
    density_dirty = false;
    density_ipack = false;
    density_norb = 0;
    density_count = 0;
    density_rho2_host = nullptr;
    density_rho2s_host = nullptr;
    density_rho2a_host = nullptr;
    density_s2_host = nullptr;
  }

  std::int64_t fallback_status() noexcept
  {
    if (!flush_density()) {
      return -1;
    }
    clear_density_state();
    return 0;
  }

  bool map_session_matches(const std::int64_t *i1_host, const double *xi1_host,
                           const std::int64_t *i3_host, const double *xi3_host,
                           std::size_t i1_count, std::size_t i3_count) const noexcept
  {
    return map_session_active && map_i1_host == i1_host && map_xi1_host == xi1_host &&
           map_i3_host == i3_host && map_xi3_host == xi3_host &&
           map_i1_count == i1_count && map_i3_count == i3_count;
  }

  bool begin_map_session(const std::int64_t *i1_host, const double *xi1_host,
                         const std::int64_t *i3_host, const double *xi3_host,
                         std::size_t i1_count, std::size_t i3_count) noexcept
  {
    if (map_session_active) {
      return false;
    }
    map_session_active = true;
    map_session_resident = false;
    map_i1_host = i1_host;
    map_xi1_host = xi1_host;
    map_i3_host = i3_host;
    map_xi3_host = xi3_host;
    map_i1_count = i1_count;
    map_i3_count = i3_count;
    return true;
  }

  void invalidate_map_residency() noexcept
  {
    map_session_resident = false;
  }

  void clear_map_session() noexcept
  {
    map_session_active = false;
    map_session_resident = false;
    map_i1_host = nullptr;
    map_xi1_host = nullptr;
    map_i3_host = nullptr;
    map_xi3_host = nullptr;
    map_i1_count = 0;
    map_i3_count = 0;
  }

  bool end_map_session() noexcept
  {
    clear_map_session();
    return true;
  }

  bool supports_grid(std::size_t x_count, std::size_t y_count, std::size_t z_count) noexcept
  {
    if (!device_ready) {
      int device = 0;
      cudaDeviceProp properties{};
      if (cudaGetDevice(&device) != cudaSuccess ||
          cudaGetDeviceProperties(&properties, device) != cudaSuccess ||
          properties.maxGridSize[0] <= 0 || properties.maxGridSize[1] <= 0 ||
          properties.maxGridSize[2] <= 0) {
        return false;
      }
      max_grid_x = static_cast<std::size_t>(properties.maxGridSize[0]);
      max_grid_y = static_cast<std::size_t>(properties.maxGridSize[1]);
      max_grid_z = static_cast<std::size_t>(properties.maxGridSize[2]);
      device_ready = true;
    }
    const std::size_t max_unsigned = static_cast<std::size_t>((std::numeric_limits<unsigned int>::max)());
    return x_count > 0 && y_count > 0 && z_count > 0 && x_count <= max_grid_x &&
           y_count <= max_grid_y && z_count <= max_grid_z && x_count <= max_unsigned &&
           y_count <= max_unsigned && z_count <= max_unsigned;
  }

  bool limit_split_count(std::size_t plane_count, std::size_t requested,
                         std::size_t *result) noexcept
  {
    if (result == nullptr || plane_count == 0 || requested == 0 ||
        !supports_grid(1, 1, 1) || plane_count > max_grid_z) {
      return false;
    }
    const std::size_t available = max_grid_z / plane_count;
    *result = requested < available ? requested : available;
    return *result > 0;
  }

  bool release() noexcept
  {
    bool success = end_session();
    success = end_map_session() && success;
    if (!x.release()) success = false;
    if (!sb.release()) success = false;
    if (!cb.release()) success = false;
    if (!xi1.release()) success = false;
    if (!xi3.release()) success = false;
    if (!xi4.release()) success = false;
    if (!xi2.release()) success = false;
    if (!i1.release()) success = false;
    if (!i3.release()) success = false;
    if (!i4.release()) success = false;
    if (!i2.release()) success = false;
    if (!rho2.release()) success = false;
    if (!rho2s.release()) success = false;
    if (!rho2a.release()) success = false;
    if (!s2_delta.release()) success = false;
    if (!staging.release()) success = false;
    device_ready = false;
    max_grid_x = 0;
    max_grid_y = 0;
    max_grid_z = 0;
    return success;
  }
};

__device__ double atomic_add(double *address, double value)
{
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

__device__ std::size_t tri_index(std::size_t a, std::size_t b)
{
  const std::size_t high = a > b ? a : b;
  const std::size_t low = a > b ? b : a;
  return high * (high - 1) / 2 + low;
}

__device__ std::size_t packed_pair_index(std::size_t first, std::size_t second)
{
  return first * (first - 1) / 2 + second;
}

__global__ void gsbbd2b_density_kernel(
    const double *x, double *rho2, double *rho2s, double *rho2a, double *s2_delta,
    std::size_t ni, std::size_t nj, std::size_t nk, std::size_t nl,
    std::size_t ioff, std::size_t joff, std::size_t koff, std::size_t loff,
    std::size_t norb, bool ipack, bool s2_active)
{
  const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  const std::size_t total = ni * nj * nk * nl;
  if (index >= total) {
    return;
  }

  std::size_t remaining = index;
  const std::size_t i = remaining % ni;
  remaining /= ni;
  const std::size_t j = remaining % nj;
  remaining /= nj;
  const std::size_t k = remaining % nk;
  const std::size_t l = remaining / nk;
  const std::size_t I = i + ioff;
  const std::size_t J = j + joff;
  const std::size_t K = k + koff;
  const std::size_t L = l + loff;
  const std::size_t ij = (J - 1) * norb + I;
  const std::size_t kl = (L - 1) * norb + K;
  const double term = (ij == kl ? 2.0 : 1.0) * x[index];

  if (!ipack) {
    atomic_add(rho2 + tri_index(ij, kl) - 1, term);
  } else {
    std::size_t i_pack = I;
    std::size_t j_pack = J;
    std::size_t k_pack = K;
    std::size_t l_pack = L;
    std::size_t ij_pack = packed_pair_index(i_pack, j_pack);
    std::size_t ji_pack = packed_pair_index(j_pack, i_pack);
    std::size_t kl_pack = packed_pair_index(k_pack, l_pack);
    double factor_pack = k_pack == l_pack ? 0.25 : 0.5;
    if (i_pack == k_pack && j_pack == l_pack) {
      factor_pack *= 0.5;
    }
    if (i_pack >= j_pack && k_pack >= l_pack && ij_pack >= kl_pack) {
      const std::size_t target = packed_pair_index(ij_pack, kl_pack) - 1;
      atomic_add(rho2s + target, factor_pack * term);
      atomic_add(rho2a + target, factor_pack * term);
    }
    if (j_pack >= i_pack && k_pack >= l_pack && ji_pack >= kl_pack) {
      const std::size_t target = packed_pair_index(ji_pack, kl_pack) - 1;
      atomic_add(rho2s + target, factor_pack * term);
      atomic_add(rho2a + target, -factor_pack * term);
    }
    i_pack = K;
    j_pack = L;
    k_pack = I;
    l_pack = J;
    ij_pack = packed_pair_index(i_pack, j_pack);
    ji_pack = packed_pair_index(j_pack, i_pack);
    kl_pack = packed_pair_index(k_pack, l_pack);
    factor_pack = k_pack == l_pack ? 0.25 : 0.5;
    if (i_pack == k_pack && j_pack == l_pack) {
      factor_pack *= 0.5;
    }
    if (i_pack >= j_pack && k_pack >= l_pack && ij_pack >= kl_pack) {
      const std::size_t target = packed_pair_index(ij_pack, kl_pack) - 1;
      atomic_add(rho2s + target, factor_pack * term);
      atomic_add(rho2a + target, factor_pack * term);
    }
    if (j_pack >= i_pack && k_pack >= l_pack && ji_pack >= kl_pack) {
      const std::size_t target = packed_pair_index(ji_pack, kl_pack) - 1;
      atomic_add(rho2s + target, factor_pack * term);
      atomic_add(rho2a + target, -factor_pack * term);
    }
  }
  if (s2_active && k == j && l == i) {
    atomic_add(s2_delta, -x[index]);
  }
}

__global__ void gsbbd2b_cuda_kernel(
    double *x, const double *sb, const double *cb, const std::int64_t *i1,
    const double *xi1, const std::int64_t *i3, const double *xi3,
    const std::int64_t *i4, const double *xi4, const std::int64_t *i2,
    const double *xi2, std::size_t nib, std::size_t njb, std::size_t nkastr,
    std::size_t kabot0, std::size_t lkabtc, std::size_t nkbstr, std::size_t ni,
    std::size_t nj, std::size_t nk, std::size_t split_count, std::size_t maps_in_split)
{
  __shared__ double a_tile[tile_size][tile_size + 1];
  __shared__ double b_tile[tile_size][tile_size + 1];
  __shared__ std::int64_t ib;
  __shared__ std::int64_t jb;
  __shared__ double factor;

  const std::size_t plane_split = static_cast<std::size_t>(blockIdx.z);
  const std::size_t plane = plane_split / split_count;
  const std::size_t split = plane_split - plane * split_count;
  const std::size_t k = plane % nk;
  const std::size_t l = plane / nk;
  const std::size_t i = static_cast<std::size_t>(blockIdx.x) * tile_size + threadIdx.y;
  const std::size_t j = static_cast<std::size_t>(blockIdx.y) * tile_size + threadIdx.x;
  const std::size_t map_begin = split * maps_in_split;
  const std::size_t map_end = map_begin + maps_in_split < nkbstr ?
                              map_begin + maps_in_split : nkbstr;
  double result = 0.0;
  for (std::size_t kb = map_begin; kb < map_end; ++kb) {
    const std::size_t k_map = kb + k * nkbstr;
    const std::size_t l_map = kb + l * nkbstr;
    if (threadIdx.x == 0 && threadIdx.y == 0) {
      ib = i4[k_map];
      jb = i2[l_map];
      factor = ib == 0 || jb == 0 ? 0.0 : xi4[k_map] * xi2[l_map];
    }
    __syncthreads();
    if (ib != 0 && jb != 0) {
      double product = 0.0;
      for (std::size_t base = 0; base < lkabtc; base += tile_size) {
        const std::size_t a_ka = base + threadIdx.x;
        const std::size_t b_ka = base + threadIdx.y;
        double a_value = 0.0;
        if (i < ni && a_ka < lkabtc) {
          const std::size_t map = kabot0 + a_ka + i * nkastr;
          const std::int64_t ia = i3[map];
          if (ia != 0) {
            a_value = xi3[map] * sb[static_cast<std::size_t>(ib - 1) +
                                     static_cast<std::size_t>(ia - 1) * nib];
          }
        }
        a_tile[threadIdx.y][threadIdx.x] = a_value;

        double b_value = 0.0;
        if (j < nj && b_ka < lkabtc) {
          const std::size_t map = kabot0 + b_ka + j * nkastr;
          const std::int64_t ja = i1[map];
          if (ja != 0) {
            b_value = xi1[map] * cb[static_cast<std::size_t>(jb - 1) +
                                     static_cast<std::size_t>(ja - 1) * njb];
          }
        }
        b_tile[threadIdx.y][threadIdx.x] = b_value;
        __syncthreads();
        if (i < ni && j < nj) {
          for (std::size_t ka = 0; ka < tile_size; ++ka) {
            product += a_tile[threadIdx.y][ka] * b_tile[ka][threadIdx.x];
          }
        }
        __syncthreads();
      }
      result += factor * product;
    }
    __syncthreads();
  }

  if (i < ni && j < nj && result != 0.0) {
    const std::size_t output = i + j * ni + (k + l * nk) * ni * nj;
    atomic_add(x + output, result);
  }
}

Workspace workspace;

} // namespace

extern "C" int64_t lucia_gsbbd2b_cuda_begin(
    const double *sb, const double *cb, double *rho2, double *rho2s, double *rho2a,
    double *s2_term1, int64_t nia, int64_t nib, int64_t nja, int64_t njb,
    int64_t norb, int64_t ipack)
{
  if (sb == nullptr || cb == nullptr || s2_term1 == nullptr || workspace.session_active ||
      (ipack != 0 && ipack != 1) || (ipack == 0 ? rho2 == nullptr : rho2s == nullptr || rho2a == nullptr)) {
    return 0;
  }
  std::size_t nia_size = 0;
  std::size_t nib_size = 0;
  std::size_t nja_size = 0;
  std::size_t njb_size = 0;
  std::size_t sb_count = 0;
  std::size_t cb_count = 0;
  std::size_t norb_size = 0;
  std::size_t norb_squared = 0;
  std::size_t density_count = 0;
  if (!lucia_cuda::positive_size(nia, &nia_size) ||
      !lucia_cuda::positive_size(nib, &nib_size) ||
      !lucia_cuda::positive_size(nja, &nja_size) ||
      !lucia_cuda::positive_size(njb, &njb_size) ||
      !lucia_cuda::positive_size(norb, &norb_size) ||
      !lucia_cuda::checked_mul(nia_size, nib_size, &sb_count) ||
      !lucia_cuda::checked_mul(nja_size, njb_size, &cb_count) ||
      !lucia_cuda::checked_mul(norb_size, norb_size, &norb_squared) ||
      (ipack == 0 ? !checked_ntri(norb_squared, &density_count)
                  : !checked_ntri(norb_size, &norb_squared) ||
                        !checked_ntri(norb_squared, &density_count)) ||
      !workspace.begin_session(sb, cb, sb_count, cb_count) ||
      !workspace.begin_density_session(rho2, rho2s, rho2a, s2_term1, norb_size,
                                       ipack == 1, density_count)) {
    workspace.clear_session();
    workspace.clear_density_state();
    return 0;
  }
  return 1;
}

extern "C" int64_t lucia_gsbbd2b_cuda_begin_maps(
    const std::int64_t *i1, const double *xi1, const std::int64_t *i3, const double *xi3,
    int64_t nkastr, int64_t ni, int64_t nj)
{
  if (i1 == nullptr || xi1 == nullptr || i3 == nullptr || xi3 == nullptr ||
      workspace.map_session_active) {
    return 0;
  }
  std::size_t nkastr_size = 0;
  std::size_t ni_size = 0;
  std::size_t nj_size = 0;
  std::size_t i1_count = 0;
  std::size_t i3_count = 0;
  if (!lucia_cuda::positive_size(nkastr, &nkastr_size) ||
      !lucia_cuda::positive_size(ni, &ni_size) ||
      !lucia_cuda::positive_size(nj, &nj_size) ||
      !lucia_cuda::checked_mul(nkastr_size, nj_size, &i1_count) ||
      !lucia_cuda::checked_mul(nkastr_size, ni_size, &i3_count) ||
      !workspace.begin_map_session(i1, xi1, i3, xi3, i1_count, i3_count)) {
    return 0;
  }
  return 1;
}

extern "C" int64_t lucia_gsbbd2b_cuda_end_maps()
{
  return workspace.end_map_session() ? 1 : -1;
}

extern "C" int64_t lucia_gsbbd2b_cuda_end()
{
  return workspace.end_session() ? 1 : -1;
}

extern "C" int64_t lucia_gsbbd2b_cuda_route(
    double *x, const double *sb, const double *cb, const int64_t *i1,
    const double *xi1, const int64_t *i3, const double *xi3, const int64_t *i4,
    const double *xi4, const int64_t *i2, const double *xi2, int64_t nia,
    int64_t nib, int64_t nja, int64_t njb, int64_t nkastr, int64_t kabot,
    int64_t lkabtc, int64_t nkbstr, int64_t ni, int64_t nj, int64_t nk,
    int64_t nl, int64_t ikord, int64_t ioff, int64_t joff, int64_t koff,
    int64_t loff, int64_t norb, int64_t ipack, int64_t s2_active)
{
  if (x == nullptr || sb == nullptr || cb == nullptr || i1 == nullptr || xi1 == nullptr ||
      i3 == nullptr || xi3 == nullptr || i4 == nullptr || xi4 == nullptr || i2 == nullptr ||
      xi2 == nullptr || ikord != 0) {
    return workspace.fallback_status();
  }
  if (workspace.session_active && !workspace.density_active) {
    return workspace.fallback_status();
  }

  std::size_t nia_size = 0, nib_size = 0, nja_size = 0, njb_size = 0;
  std::size_t nkastr_size = 0, kabot_size = 0, lkabtc_size = 0, nkbstr_size = 0;
  std::size_t ni_size = 0, nj_size = 0, nk_size = 0, nl_size = 0;
  std::size_t ioff_size = 0, joff_size = 0, koff_size = 0, loff_size = 0;
  std::size_t norb_size = 0;
  if (!lucia_cuda::positive_size(nia, &nia_size) || !lucia_cuda::positive_size(nib, &nib_size) ||
      !lucia_cuda::positive_size(nja, &nja_size) || !lucia_cuda::positive_size(njb, &njb_size) ||
      !lucia_cuda::positive_size(nkastr, &nkastr_size) || !lucia_cuda::positive_size(kabot, &kabot_size) ||
      !lucia_cuda::positive_size(lkabtc, &lkabtc_size) || !lucia_cuda::positive_size(nkbstr, &nkbstr_size) ||
      !lucia_cuda::positive_size(ni, &ni_size) || !lucia_cuda::positive_size(nj, &nj_size) ||
      !lucia_cuda::positive_size(nk, &nk_size) || !lucia_cuda::positive_size(nl, &nl_size) ||
      !lucia_cuda::positive_size(ioff, &ioff_size) || !lucia_cuda::positive_size(joff, &joff_size) ||
      !lucia_cuda::positive_size(koff, &koff_size) || !lucia_cuda::positive_size(loff, &loff_size) ||
      !lucia_cuda::positive_size(norb, &norb_size)) {
    return workspace.fallback_status();
  }
  if ((ipack != 0 && ipack != 1) || (s2_active != 0 && s2_active != 1) ||
      !valid_extent(ioff_size, ni_size, norb_size) ||
      !valid_extent(joff_size, nj_size, norb_size) ||
      !valid_extent(koff_size, nk_size, norb_size) ||
      !valid_extent(loff_size, nl_size, norb_size) ||
      (s2_active == 1 && (nk_size != nj_size || nl_size != ni_size))) {
    return workspace.fallback_status();
  }
  if (workspace.density_active &&
      (norb_size != workspace.density_norb || (ipack == 1) != workspace.density_ipack)) {
    return workspace.fallback_status();
  }
  const std::size_t kabot0 = kabot_size - 1;
  if (kabot0 >= nkastr_size || lkabtc_size > nkastr_size - kabot0) {
    return workspace.fallback_status();
  }

  std::size_t x_count = 0, sb_count = 0, cb_count = 0;
  std::size_t i1_count = 0, i3_count = 0, i4_count = 0, i2_count = 0;
  std::size_t partial = 0;
  if (!lucia_cuda::checked_mul(ni_size, nj_size, &partial) ||
      !lucia_cuda::checked_mul(partial, nk_size, &partial) ||
      !lucia_cuda::checked_mul(partial, nl_size, &x_count) ||
      !lucia_cuda::checked_mul(nia_size, nib_size, &sb_count) ||
      !lucia_cuda::checked_mul(nja_size, njb_size, &cb_count) ||
      !lucia_cuda::checked_mul(nkastr_size, nj_size, &i1_count) ||
      !lucia_cuda::checked_mul(nkastr_size, ni_size, &i3_count) ||
      !lucia_cuda::checked_mul(nkbstr_size, nk_size, &i4_count) ||
      !lucia_cuda::checked_mul(nkbstr_size, nl_size, &i2_count)) {
    return workspace.fallback_status();
  }

  const bool session_match = workspace.session_matches(sb, cb, sb_count, cb_count);
  if (workspace.session_active && !session_match) {
    return workspace.fallback_status();
  }
  const bool shared_match = lucia_sigma_cuda_blocks::sigma_blocks_match(
      const_cast<double *>(sb), cb, sb_count, cb_count);
  const bool map_match = workspace.map_session_matches(i1, xi1, i3, xi3, i1_count, i3_count);
  if (workspace.map_session_active && !map_match) {
    return workspace.fallback_status();
  }

  std::size_t norb_squared = 0;
  std::size_t density_count = 0;
  if (!lucia_cuda::checked_mul(norb_size, norb_size, &norb_squared) ||
      (ipack == 0 ? !checked_ntri(norb_squared, &density_count)
                  : !checked_ntri(norb_size, &norb_squared) ||
                        !checked_ntri(norb_squared, &density_count))) {
    return workspace.fallback_status();
  }
  if (workspace.density_active && density_count != workspace.density_count) {
    return workspace.fallback_status();
  }

  std::size_t active_pairs = 0;
  for (std::size_t kb = 0; kb < nkbstr_size; ++kb) {
    std::size_t active_k = 0;
    for (std::size_t k = 0; k < nk_size; ++k) {
      if (i4[kb + k * nkbstr_size] != 0) ++active_k;
    }
    std::size_t active_l = 0;
    for (std::size_t l = 0; l < nl_size; ++l) {
      if (i2[kb + l * nkbstr_size] != 0) ++active_l;
    }
    std::size_t pairs = 0;
    if (!lucia_cuda::checked_mul(active_k, active_l, &pairs) ||
        !checked_add(active_pairs, pairs, &active_pairs)) {
      return workspace.fallback_status();
    }
  }
  std::size_t active_work = active_pairs;
  if (!lucia_cuda::checked_mul(active_work, ni_size, &active_work) ||
      !lucia_cuda::checked_mul(active_work, nj_size, &active_work) ||
      !lucia_cuda::checked_mul(active_work, lkabtc_size, &active_work) ||
      active_work < minimum_active_work) {
    return workspace.fallback_status();
  }

  for (std::size_t j = 0; j < nj_size; ++j) {
    for (std::size_t ka = 0; ka < lkabtc_size; ++ka) {
      const std::size_t map = kabot0 + ka + j * nkastr_size;
      if (i1[map] < 0 || i1[map] > nja || !std::isfinite(xi1[map])) {
        return workspace.fallback_status();
      }
    }
  }
  for (std::size_t i = 0; i < ni_size; ++i) {
    for (std::size_t ka = 0; ka < lkabtc_size; ++ka) {
      const std::size_t map = kabot0 + ka + i * nkastr_size;
      if (i3[map] < 0 || i3[map] > nia || !std::isfinite(xi3[map])) {
        return workspace.fallback_status();
      }
    }
  }
  for (std::size_t p = 0; p < i4_count; ++p) {
    if (i4[p] < 0 || i4[p] > nib || !std::isfinite(xi4[p])) {
      return workspace.fallback_status();
    }
  }
  for (std::size_t p = 0; p < i2_count; ++p) {
    if (i2[p] < 0 || i2[p] > njb || !std::isfinite(xi2[p])) {
      return workspace.fallback_status();
    }
  }

  const std::size_t grid_x = ni_size / tile_size + (ni_size % tile_size != 0);
  const std::size_t grid_y = nj_size / tile_size + (nj_size % tile_size != 0);
  const std::size_t requested_splits =
      nkbstr_size / maps_per_split + (nkbstr_size % maps_per_split != 0);
  std::size_t plane_count = 0, split_count = 0, grid_z = 0;
  if (!lucia_cuda::checked_mul(nk_size, nl_size, &plane_count) ||
      !workspace.limit_split_count(plane_count, requested_splits, &split_count) ||
      !lucia_cuda::checked_mul(plane_count, split_count, &grid_z) ||
      !workspace.supports_grid(grid_x, grid_y, grid_z)) {
    return workspace.fallback_status();
  }
  const std::size_t maps_in_split =
      nkbstr_size / split_count + (nkbstr_size % split_count != 0);

  constexpr std::size_t density_block_size = 256;
  std::size_t density_grid = x_count / density_block_size;
  if (x_count % density_block_size != 0 && !checked_add(density_grid, 1, &density_grid)) {
    return workspace.fallback_status();
  }
  if (workspace.density_active &&
      (!workspace.density_uploaded || !workspace.supports_grid(density_grid, 1, 1))) {
    return workspace.fallback_status();
  }

  std::size_t x_bytes = 0;
  if (!lucia_cuda::checked_mul(x_count, sizeof(double), &x_bytes) ||
      !workspace.x.reserve(x_count)) {
    return workspace.fallback_status();
  }
  if (cudaMemset(workspace.x.data, 0, x_bytes) != cudaSuccess) {
    workspace.invalidate_map_residency();
    return -1;
  }
  const bool upload_blocks = session_match && !workspace.session_resident;
  double *device_sb = workspace.sb.data;
  const double *device_cb = workspace.cb.data;
  if (shared_match && !lucia_sigma_cuda_blocks::sigma_blocks_acquire(
                          const_cast<double *>(sb), cb, sb_count, cb_count,
                          &device_sb, &device_cb)) {
    return workspace.fallback_status();
  }
  if (!shared_match && (!session_match || upload_blocks)) {
    int upload_status = upload_device(workspace.sb, sb, sb_count);
    if (upload_status == 1) {
      upload_status = upload_device(workspace.cb, cb, cb_count);
    }
    if (upload_status != 1) {
      workspace.invalidate_map_residency();
      return upload_status < 0 ? -1 : workspace.fallback_status();
    }
  }
  if (!shared_match) {
    device_sb = workspace.sb.data;
    device_cb = workspace.cb.data;
  }
  if (!shared_match && session_match && upload_blocks) {
    workspace.session_resident = true;
  }
  const bool upload_maps = map_match && !workspace.map_session_resident;
  if (!map_match || upload_maps) {
    int upload_status = upload_device(workspace.i1, i1, i1_count);
    if (upload_status == 1) {
      upload_status = upload_device(workspace.xi1, xi1, i1_count);
    }
    if (upload_status == 1) {
      upload_status = upload_device(workspace.i3, i3, i3_count);
    }
    if (upload_status == 1) {
      upload_status = upload_device(workspace.xi3, xi3, i3_count);
    }
    if (upload_status != 1) {
      workspace.invalidate_map_residency();
      return upload_status < 0 ? -1 : workspace.fallback_status();
    }
  }
  if (map_match && upload_maps) {
    workspace.map_session_resident = true;
  }
  int upload_status = upload_device(workspace.i4, i4, i4_count);
  if (upload_status == 1) {
    upload_status = upload_device(workspace.xi4, xi4, i4_count);
  }
  if (upload_status == 1) {
    upload_status = upload_device(workspace.i2, i2, i2_count);
  }
  if (upload_status == 1) {
    upload_status = upload_device(workspace.xi2, xi2, i2_count);
  }
  if (upload_status != 1) {
    workspace.invalidate_map_residency();
    return upload_status < 0 ? -1 : workspace.fallback_status();
  }
  if (!workspace.density_active && !workspace.staging.reserve(x_count)) {
    workspace.invalidate_map_residency();
    return workspace.fallback_status();
  }

  const dim3 block(static_cast<unsigned int>(tile_size), static_cast<unsigned int>(tile_size), 1);
  const dim3 grid(static_cast<unsigned int>(grid_x), static_cast<unsigned int>(grid_y),
                  static_cast<unsigned int>(grid_z));
  gsbbd2b_cuda_kernel<<<grid, block>>>(
      workspace.x.data, device_sb, device_cb, workspace.i1.data,
      workspace.xi1.data, workspace.i3.data, workspace.xi3.data, workspace.i4.data,
      workspace.xi4.data, workspace.i2.data, workspace.xi2.data, nib_size, njb_size,
      nkastr_size, kabot0, lkabtc_size, nkbstr_size, ni_size, nj_size, nk_size,
      split_count, maps_in_split);
  if (cudaGetLastError() != cudaSuccess) {
    workspace.invalidate_map_residency();
    return -1;
  }
  if (workspace.density_active) {
    gsbbd2b_density_kernel<<<dim3(static_cast<unsigned int>(density_grid)),
                             dim3(static_cast<unsigned int>(density_block_size))>>>(
        workspace.x.data, workspace.rho2.data, workspace.rho2s.data, workspace.rho2a.data,
        workspace.s2_delta.data, ni_size, nj_size, nk_size, nl_size, ioff_size, joff_size,
        koff_size, loff_size, norb_size, ipack == 1, s2_active == 1);
    if (cudaGetLastError() != cudaSuccess) {
      workspace.invalidate_map_residency();
      return -1;
    }
    workspace.density_dirty = true;
    return 1;
  }
  if (cudaDeviceSynchronize() != cudaSuccess ||
      cudaMemcpy(workspace.staging.data, workspace.x.data, x_bytes,
                 cudaMemcpyDeviceToHost) != cudaSuccess) {
    workspace.invalidate_map_residency();
    return -1;
  }
  std::memcpy(x, workspace.staging.data, x_bytes);
  return 1;
}

extern "C" void lucia_gsbbd2b_cuda_release()
{
  workspace.release();
}
