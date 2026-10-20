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

constexpr std::size_t block_size = 256;
constexpr std::size_t tile_size = 16;

bool checked_add(std::size_t left, std::size_t right, std::size_t *result) noexcept {
  if (left > (std::numeric_limits<std::size_t>::max)() - right) {
    return false;
  }
  *result = left + right;
  return true;
}
bool append_section(std::size_t bytes, std::size_t *offset, std::size_t *total) noexcept {
  *offset = *total;
  return checked_add(*total, bytes, total);
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
  lucia_cuda::DeviceBuffer<double> sb;
  lucia_cuda::DeviceBuffer<double> cb;
  lucia_cuda::DeviceBuffer<double> xint;
  lucia_cuda::DeviceBuffer<double> product;
  lucia_cuda::PinnedBuffer<double> staging;
  RouteSlot route[2];
  unsigned int next_route = 0;
  bool device_ready = false;
  bool session_active = false;
  bool session_resident = false;
  double *session_sb_host = nullptr;
  const double *session_cb_host = nullptr;
  std::size_t session_sb_count = 0;
  std::size_t session_cb_count = 0;
  bool xint_session_active = false;
  bool xint_session_resident = false;
  const double *xint_session_host = nullptr;
  std::size_t xint_session_count = 0;
  std::size_t max_grid_x = 0;
  std::size_t max_grid_y = 0;
  std::size_t max_grid_z = 0;

  bool session_matches(double *sb_host, const double *cb_host, std::size_t sb_count, std::size_t cb_count) const noexcept {
    return session_active && session_sb_host == sb_host && session_cb_host == cb_host && session_sb_count == sb_count
           && session_cb_count == cb_count;
  }

  bool begin_session(double *sb_host, const double *cb_host, std::size_t sb_count, std::size_t cb_count) noexcept {
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

  bool xint_session_matches(const double *xint_host, std::size_t xint_count) const noexcept {
    return xint_session_active && xint_session_host == xint_host && xint_session_count == xint_count;
  }

  bool begin_xint_session(const double *xint_host, std::size_t xint_count) noexcept {
    if (xint_session_active) {
      return false;
    }
    xint_session_active = true;
    xint_session_resident = false;
    xint_session_host = xint_host;
    xint_session_count = xint_count;
    return true;
  }

  void invalidate_xint_residency() noexcept {
    xint_session_resident = false;
  }

  void clear_xint_session() noexcept {
    xint_session_active = false;
    xint_session_resident = false;
    xint_session_host = nullptr;
    xint_session_count = 0;
  }

  bool end_xint_session() noexcept {
    clear_xint_session();
    return true;
  }

  bool flush_session() noexcept {
    if (!session_active || !session_resident) {
      return true;
    }
    invalidate_residency();
    std::size_t sb_bytes = 0;
    if (!lucia_cuda::checked_mul(session_sb_count, sizeof(double), &sb_bytes) || !staging.reserve(session_sb_count)
        || cudaMemcpy(staging.data, sb.data, sb_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) {
      return false;
    }
    std::memcpy(session_sb_host, staging.data, sb_bytes);
    return true;
  }

  void invalidate_residency() noexcept {
    session_resident = false;
  }

  std::int64_t fallback_status() noexcept {
    return flush_session() ? 0 : -1;
  }

  void clear_session() noexcept {
    session_active = false;
    session_resident = false;
    session_sb_host = nullptr;
    session_cb_host = nullptr;
    session_sb_count = 0;
    session_cb_count = 0;
  }

  bool end_session() noexcept {
    const bool success = flush_session();
    clear_session();
    return success;
  }

  bool supports_grid(std::size_t x, std::size_t y, std::size_t z) noexcept {
    if (!device_ready) {
      int device = 0;
      cudaDeviceProp properties{};
      if (cudaGetDevice(&device) != cudaSuccess || cudaGetDeviceProperties(&properties, device) != cudaSuccess) {
        return false;
      }
      max_grid_x = static_cast<std::size_t>(properties.maxGridSize[0]);
      max_grid_y = static_cast<std::size_t>(properties.maxGridSize[1]);
      max_grid_z = static_cast<std::size_t>(properties.maxGridSize[2]);
      device_ready = true;
    }
    return x > 0 && y > 0 && z > 0 && x <= max_grid_x && y <= max_grid_y && z <= max_grid_z;
  }

  bool release() noexcept {
    bool success = end_session();
    success = end_xint_session() && success;
    if (!sb.release())
      success = false;
    if (!cb.release())
      success = false;
    if (!xint.release())
      success = false;
    if (!product.release())
      success = false;
    if (!staging.release())
      success = false;
    if (!route[0].release())
      success = false;
    if (!route[1].release())
      success = false;
    next_route = 0;
    device_ready = false;
    max_grid_x = 0;
    max_grid_y = 0;
    max_grid_z = 0;
    return success;
  }
};

std::int64_t route_fallback(Workspace &workspace, double *sb, const double *cb) noexcept {
  bool matched = false;
  if (!lucia_sigma_cuda_blocks::sigma_blocks_flush_if_match(sb, cb, &matched)) {
    return -1;
  }
  return matched ? 0 : workspace.fallback_status();
}

__device__ double atomic_add(double *address, double value);

__global__ void rsbb2a_contract_kernel(double *product, const double *cb, const double *xint, const std::int64_t *active_source,
                                       const std::int32_t *active_jl, const std::int32_t *active_offset, const double *active_sign,
                                       std::size_t nrow, std::size_t ibot0, std::size_t nibtc, std::size_t nkbtc, std::size_t nik,
                                       double factor) {
  __shared__ double a[tile_size][tile_size + 1];
  __shared__ double b[tile_size][tile_size + 1];

  const std::size_t row = static_cast<std::size_t>(blockIdx.x) * tile_size + threadIdx.x;
  const std::size_t ik = static_cast<std::size_t>(blockIdx.y) * tile_size + threadIdx.y;
  const std::size_t k = blockIdx.z;
  const std::int32_t begin = active_offset[k];
  const std::int32_t count = active_offset[k + 1] - begin;
  double sum = 0.0;

  for (std::int32_t base = 0; base < count; base += tile_size) {
    const std::int32_t entry = begin + base + threadIdx.y;
    double a_value = 0.0;
    if (row < nibtc && base + threadIdx.y < count) {
      a_value = active_sign[entry] * cb[ibot0 + row + static_cast<std::size_t>(active_source[entry] - 1) * nrow];
    }
    a[threadIdx.x][threadIdx.y] = a_value;

    const std::size_t b_ik = static_cast<std::size_t>(blockIdx.y) * tile_size + threadIdx.x;
    const std::int32_t b_entry = begin + base + threadIdx.y;
    b[threadIdx.x][threadIdx.y]
        = b_ik < nik && base + threadIdx.y < count ? xint[b_ik + static_cast<std::size_t>(active_jl[b_entry]) * nik] : 0.0;
    __syncthreads();

    if (row < nibtc && ik < nik) {
#pragma unroll
      for (std::size_t q = 0; q < tile_size; ++q) {
        sum += a[threadIdx.x][q] * b[threadIdx.y][q];
      }
    }
    __syncthreads();
  }

  if (row < nibtc && ik < nik) {
    const std::size_t plane = nibtc * nkbtc;
    product[row + k * nibtc + ik * plane] = factor * sum;
  }
}

__global__ void rsbb2a_scatter_kernel(double *sb, const double *product, const std::int64_t *smap, const double *ssign,
                                      std::size_t nrow, std::size_t ibot0, std::size_t nibtc, std::size_t nkbtc, std::size_t nik) {
  const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  const std::size_t plane = nibtc * nkbtc;
  if (index >= plane * nik) {
    return;
  }

  const std::size_t ik = index / plane;
  const std::size_t within = index - ik * plane;
  const std::size_t k = within / nibtc;
  const std::size_t row = within - k * nibtc;
  const std::size_t map_offset = k + ik * nkbtc;
  const std::int64_t destination = smap[map_offset];
  if (destination != 0) {
    const std::size_t sb_offset = ibot0 + row + static_cast<std::size_t>(destination - 1) * nrow;
    atomic_add(sb + sb_offset, ssign[map_offset] * product[index]);
  }
}

Workspace workspace;

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

} // namespace

extern "C" int64_t lucia_rsbb2a_cuda_begin(double *sb, const double *cb, int64_t nrow, int64_t nsb, int64_t ncb) {
  if (sb == nullptr || cb == nullptr) {
    return 0;
  }
  std::size_t nrow_size = 0;
  std::size_t nsb_size = 0;
  std::size_t ncb_size = 0;
  std::size_t sb_count = 0;
  std::size_t cb_count = 0;
  if (!lucia_cuda::positive_size(nrow, &nrow_size) || !lucia_cuda::positive_size(nsb, &nsb_size)
      || !lucia_cuda::positive_size(ncb, &ncb_size) || !lucia_cuda::checked_mul(nrow_size, nsb_size, &sb_count)
      || !lucia_cuda::checked_mul(nrow_size, ncb_size, &cb_count)) {
    return 0;
  }
  if (lucia_sigma_cuda_blocks::sigma_blocks_match(sb, cb, sb_count, cb_count)) {
    return 1;
  }
  if (workspace.session_active || !workspace.begin_session(sb, cb, sb_count, cb_count)) {
    return 0;
  }
  return 1;
}

extern "C" int64_t lucia_rsbb2a_cuda_xint_begin(const double *xint, int64_t nik, int64_t njl) {
  if (xint == nullptr) {
    return 0;
  }
  std::size_t nik_size = 0;
  std::size_t njl_size = 0;
  std::size_t xint_count = 0;
  if (!lucia_cuda::positive_size(nik, &nik_size) || !lucia_cuda::positive_size(njl, &njl_size)
      || !lucia_cuda::checked_mul(nik_size, njl_size, &xint_count) || !workspace.begin_xint_session(xint, xint_count)) {
    return 0;
  }
  return 1;
}

extern "C" int64_t lucia_rsbb2a_cuda_xint_end() {
  return workspace.end_xint_session() ? 1 : -1;
}

extern "C" int64_t lucia_rsbb2a_cuda_end() {
  return workspace.end_session() ? 1 : -1;
}

extern "C" int64_t lucia_rsbb2a_cuda_route(double *sb, const double *cb, const double *xint, const int64_t *cmap,
                                           const double *csign, const int64_t *smap, const double *ssign, int64_t nrow, int64_t nsb,
                                           int64_t ncb, int64_t ibot, int64_t nibtc, int64_t nkbtc, int64_t nik, int64_t njl,
                                           double factor) {
  if (sb == nullptr || cb == nullptr || xint == nullptr || cmap == nullptr || csign == nullptr || smap == nullptr
      || ssign == nullptr || !std::isfinite(factor)) {
    return route_fallback(workspace, sb, cb);
  }

  std::size_t nrow_size = 0;
  std::size_t nsb_size = 0;
  std::size_t ncb_size = 0;
  std::size_t ibot_size = 0;
  std::size_t nibtc_size = 0;
  std::size_t nkbtc_size = 0;
  std::size_t nik_size = 0;
  std::size_t njl_size = 0;
  if (!lucia_cuda::positive_size(nrow, &nrow_size) || !lucia_cuda::positive_size(nsb, &nsb_size)
      || !lucia_cuda::positive_size(ncb, &ncb_size) || !lucia_cuda::positive_size(ibot, &ibot_size)
      || !lucia_cuda::positive_size(nibtc, &nibtc_size) || !lucia_cuda::positive_size(nkbtc, &nkbtc_size)
      || !lucia_cuda::positive_size(nik, &nik_size) || !lucia_cuda::positive_size(njl, &njl_size)) {
    return route_fallback(workspace, sb, cb);
  }

  const std::size_t ibot0 = ibot_size - 1;
  if (ibot0 >= nrow_size || nibtc_size > nrow_size - ibot0) {
    return route_fallback(workspace, sb, cb);
  }

  std::size_t sb_count = 0;
  std::size_t cb_count = 0;
  std::size_t c_map_count = 0;
  std::size_t s_map_count = 0;
  std::size_t xint_count = 0;
  std::size_t plane = 0;
  std::size_t product_count = 0;
  if (!lucia_cuda::checked_mul(nrow_size, nsb_size, &sb_count) || !lucia_cuda::checked_mul(nrow_size, ncb_size, &cb_count)
      || !lucia_cuda::checked_mul(nkbtc_size, njl_size, &c_map_count)
      || !lucia_cuda::checked_mul(nkbtc_size, nik_size, &s_map_count) || !lucia_cuda::checked_mul(nik_size, njl_size, &xint_count)
      || !lucia_cuda::checked_mul(nibtc_size, nkbtc_size, &plane) || !lucia_cuda::checked_mul(plane, nik_size, &product_count)) {
    return route_fallback(workspace, sb, cb);
  }

  const bool xint_match = workspace.xint_session_matches(xint, xint_count);
  if (workspace.xint_session_active && !xint_match) {
    return route_fallback(workspace, sb, cb);
  }
  const bool shared_match = lucia_sigma_cuda_blocks::sigma_blocks_match(sb, cb, sb_count, cb_count);
  const bool local_match = workspace.session_matches(sb, cb, sb_count, cb_count);
  if (workspace.session_active && !local_match && !shared_match && route_fallback(workspace, sb, cb) != 0) {
    return -1;
  }

  std::size_t active_count = 0;
  for (std::size_t i = 0; i < c_map_count; ++i) {
    if (cmap[i] < 0 || cmap[i] > ncb || !std::isfinite(csign[i])) {
      return route_fallback(workspace, sb, cb);
    }
    if (cmap[i] != 0) {
      ++active_count;
    }
  }
  for (std::size_t i = 0; i < s_map_count; ++i) {
    if (smap[i] < 0 || smap[i] > nsb || !std::isfinite(ssign[i])) {
      return route_fallback(workspace, sb, cb);
    }
  }

  std::size_t offset_count = 0;
  const std::size_t max_compact = static_cast<std::size_t>((std::numeric_limits<std::int32_t>::max)());
  std::size_t source_bytes = 0;
  std::size_t sign_bytes = 0;
  std::size_t smap_bytes = 0;
  std::size_t ssign_bytes = 0;
  std::size_t jl_bytes = 0;
  std::size_t offset_bytes = 0;
  std::size_t source_offset = 0;
  std::size_t sign_offset = 0;
  std::size_t smap_offset = 0;
  std::size_t ssign_offset = 0;
  std::size_t jl_offset = 0;
  std::size_t offset_offset = 0;
  std::size_t payload_bytes = 0;
  if (c_map_count > max_compact || njl_size > max_compact || !checked_add(nkbtc_size, 1, &offset_count)
      || !lucia_cuda::checked_mul(active_count, sizeof(std::int64_t), &source_bytes)
      || !lucia_cuda::checked_mul(active_count, sizeof(double), &sign_bytes)
      || !lucia_cuda::checked_mul(s_map_count, sizeof(std::int64_t), &smap_bytes)
      || !lucia_cuda::checked_mul(s_map_count, sizeof(double), &ssign_bytes)
      || !lucia_cuda::checked_mul(active_count, sizeof(std::int32_t), &jl_bytes)
      || !lucia_cuda::checked_mul(offset_count, sizeof(std::int32_t), &offset_bytes)
      || !append_section(source_bytes, &source_offset, &payload_bytes) || !append_section(sign_bytes, &sign_offset, &payload_bytes)
      || !append_section(smap_bytes, &smap_offset, &payload_bytes) || !append_section(ssign_bytes, &ssign_offset, &payload_bytes)
      || !append_section(jl_bytes, &jl_offset, &payload_bytes) || !append_section(offset_bytes, &offset_offset, &payload_bytes)) {
    return route_fallback(workspace, sb, cb);
  }

  RouteSlot &slot = workspace.route[workspace.next_route];
  if (!slot.prepare() || !slot.host.reserve(payload_bytes)) {
    return route_fallback(workspace, sb, cb);
  }
  unsigned char *payload_host = slot.host.data;
  std::int64_t *active_source_host = reinterpret_cast<std::int64_t *>(payload_host + source_offset);
  double *active_sign_host = reinterpret_cast<double *>(payload_host + sign_offset);
  std::int64_t *smap_host = reinterpret_cast<std::int64_t *>(payload_host + smap_offset);
  double *ssign_host = reinterpret_cast<double *>(payload_host + ssign_offset);
  std::int32_t *active_jl_host = reinterpret_cast<std::int32_t *>(payload_host + jl_offset);
  std::int32_t *active_offset_host = reinterpret_cast<std::int32_t *>(payload_host + offset_offset);
  std::memcpy(smap_host, smap, smap_bytes);
  std::memcpy(ssign_host, ssign, ssign_bytes);
  std::size_t active_index = 0;
  for (std::size_t k = 0; k < nkbtc_size; ++k) {
    active_offset_host[k] = static_cast<std::int32_t>(active_index);
    for (std::size_t jl = 0; jl < njl_size; ++jl) {
      const std::size_t map = k + jl * nkbtc_size;
      if (cmap[map] != 0) {
        active_source_host[active_index] = cmap[map];
        active_jl_host[active_index] = static_cast<std::int32_t>(jl);
        active_sign_host[active_index] = csign[map];
        ++active_index;
      }
    }
  }
  active_offset_host[nkbtc_size] = static_cast<std::int32_t>(active_index);

  const std::size_t contract_grid_x = nibtc_size / tile_size + (nibtc_size % tile_size != 0);
  const std::size_t contract_grid_y = nik_size / tile_size + (nik_size % tile_size != 0);
  std::size_t scatter_blocks = product_count / block_size;
  if (product_count % block_size != 0 && !checked_add(scatter_blocks, 1, &scatter_blocks)) {
    return route_fallback(workspace, sb, cb);
  }
  if (!workspace.supports_grid(contract_grid_x, contract_grid_y, nkbtc_size) || !workspace.supports_grid(scatter_blocks, 1, 1)) {
    return route_fallback(workspace, sb, cb);
  }

  double *device_sb = workspace.sb.data;
  const double *device_cb = workspace.cb.data;
  if (shared_match) {
    if (!lucia_sigma_cuda_blocks::sigma_blocks_acquire(sb, cb, sb_count, cb_count, &device_sb, &device_cb)) {
      return route_fallback(workspace, sb, cb);
    }
  }
  const bool upload_blocks = local_match && !workspace.session_resident;
  if (!shared_match && (!local_match || upload_blocks)
      && (!workspace.sb.copy_from(sb, sb_count) || !workspace.cb.copy_from(cb, cb_count))) {
    return route_fallback(workspace, sb, cb);
  }
  if (!shared_match) {
    device_sb = workspace.sb.data;
    device_cb = workspace.cb.data;
  }
  const bool upload_xint = xint_match && !workspace.xint_session_resident;
  if ((!xint_match || upload_xint) && !workspace.xint.copy_from(xint, xint_count)) {
    return route_fallback(workspace, sb, cb);
  }
  if (xint_match && upload_xint) {
    workspace.xint_session_resident = true;
  }
  if (!workspace.product.reserve(product_count) || (!shared_match && !local_match && !workspace.staging.reserve(sb_count))) {
    return route_fallback(workspace, sb, cb);
  }
  if (!slot.upload(payload_bytes)) {
    return route_fallback(workspace, sb, cb);
  }
  workspace.next_route = (workspace.next_route + 1U) % 2U;

  const unsigned char *payload_device = slot.device.data;
  const std::int64_t *active_source = reinterpret_cast<const std::int64_t *>(payload_device + source_offset);
  const double *active_sign = reinterpret_cast<const double *>(payload_device + sign_offset);
  const std::int64_t *smap_device = reinterpret_cast<const std::int64_t *>(payload_device + smap_offset);
  const double *ssign_device = reinterpret_cast<const double *>(payload_device + ssign_offset);
  const std::int32_t *active_jl = reinterpret_cast<const std::int32_t *>(payload_device + jl_offset);
  const std::int32_t *active_offset = reinterpret_cast<const std::int32_t *>(payload_device + offset_offset);

  const dim3 contract_block(static_cast<unsigned int>(tile_size), static_cast<unsigned int>(tile_size));
  const dim3 contract_grid(static_cast<unsigned int>(contract_grid_x), static_cast<unsigned int>(contract_grid_y),
                           static_cast<unsigned int>(nkbtc_size));
  rsbb2a_contract_kernel<<<contract_grid, contract_block>>>(workspace.product.data, device_cb, workspace.xint.data, active_source,
                                                            active_jl, active_offset, active_sign, nrow_size, ibot0, nibtc_size,
                                                            nkbtc_size, nik_size, factor);
  bool success = cudaGetLastError() == cudaSuccess;
  if (success) {
    rsbb2a_scatter_kernel<<<static_cast<unsigned int>(scatter_blocks), block_size>>>(
        device_sb, workspace.product.data, smap_device, ssign_device, nrow_size, ibot0, nibtc_size, nkbtc_size, nik_size);
    success = cudaGetLastError() == cudaSuccess;
  }
  if (!shared_match && !local_match && cudaDeviceSynchronize() != cudaSuccess) {
    success = false;
  }
  if (!success) {
    workspace.invalidate_residency();
    workspace.invalidate_xint_residency();
    return -1;
  }
  if (shared_match) {
    return 1;
  }
  if (local_match) {
    workspace.session_resident = true;
    return 1;
  }

  std::size_t sb_bytes = 0;
  if (!lucia_cuda::checked_mul(sb_count, sizeof(double), &sb_bytes)
      || cudaMemcpy(workspace.staging.data, device_sb, sb_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) {
    return -1;
  }
  std::memcpy(sb, workspace.staging.data, sb_bytes);
  return 1;
}

extern "C" void lucia_rsbb2a_cuda_release() {
  workspace.release();
}
