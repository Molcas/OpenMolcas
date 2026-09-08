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

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>

namespace {

constexpr std::size_t threads_per_block = 128;
constexpr std::size_t minimum_traversal_work = 5000000;
constexpr std::size_t minimum_work_per_transfer_byte = 4;

static_assert(sizeof(double) == sizeof(std::int64_t), "RSBB1E transfer payload types must have equal size");

bool checked_add(std::size_t left, std::size_t right, std::size_t *result) noexcept {
  if (left > (std::numeric_limits<std::size_t>::max)() - right) {
    return false;
  }
  *result = left + right;
  return true;
}
struct Workspace {
  lucia_cuda::DeviceBuffer<double> sb;
  lucia_cuda::DeviceBuffer<double> cb;
  lucia_cuda::DeviceBuffer<double> h;
  lucia_cuda::DeviceBuffer<std::int64_t> i1;
  lucia_cuda::DeviceBuffer<double> xi1;
  lucia_cuda::DeviceBuffer<std::int64_t> i2;
  lucia_cuda::DeviceBuffer<double> xi2;
  lucia_cuda::PinnedBuffer<double> staging;
  bool device_ready = false;
  bool session_active = false;
  bool session_resident = false;
  double *session_sb_host = nullptr;
  const double *session_cb_host = nullptr;
  std::size_t session_nrow = 0;
  std::size_t session_sb_count = 0;
  std::size_t session_cb_count = 0;
  std::size_t max_grid_x = 0;

  bool session_matches(double *sb_host, const double *cb_host, std::size_t nrow) const noexcept {
    return session_active && session_sb_host == sb_host && session_cb_host == cb_host && session_nrow == nrow;
  }

  bool begin_session(double *sb_host, const double *cb_host, std::size_t nrow) noexcept {
    if (session_active) {
      return false;
    }
    session_active = true;
    session_resident = false;
    session_sb_host = sb_host;
    session_cb_host = cb_host;
    session_nrow = nrow;
    session_sb_count = 0;
    session_cb_count = 0;
    return true;
  }

  void invalidate_residency() noexcept {
    session_resident = false;
    session_sb_count = 0;
    session_cb_count = 0;
  }

  bool flush_session() noexcept {
    if (!session_active || !session_resident) {
      return true;
    }
    std::size_t sb_bytes = 0;
    if (!lucia_cuda::checked_mul(session_sb_count, sizeof(double), &sb_bytes) || !staging.reserve(session_sb_count)
        || cudaMemcpy(staging.data, sb.data, sb_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) {
      invalidate_residency();
      return false;
    }
    std::memcpy(session_sb_host, staging.data, sb_bytes);
    invalidate_residency();
    return true;
  }

  std::int64_t fallback_status() noexcept {
    return flush_session() ? 0 : -1;
  }

  void clear_session() noexcept {
    session_active = false;
    session_resident = false;
    session_sb_host = nullptr;
    session_cb_host = nullptr;
    session_nrow = 0;
    session_sb_count = 0;
    session_cb_count = 0;
  }

  bool end_session() noexcept {
    const bool success = flush_session();
    clear_session();
    return success;
  }

  bool supports_grid(std::size_t count) noexcept {
    if (!device_ready) {
      int device = 0;
      cudaDeviceProp properties{};
      if (cudaGetDevice(&device) != cudaSuccess || cudaGetDeviceProperties(&properties, device) != cudaSuccess
          || properties.maxGridSize[0] <= 0) {
        return false;
      }
      max_grid_x = static_cast<std::size_t>(properties.maxGridSize[0]);
      device_ready = true;
    }
    return count <= max_grid_x && count <= static_cast<std::size_t>((std::numeric_limits<unsigned int>::max)());
  }

  bool release() noexcept {
    bool success = end_session();
    if (!sb.release())
      success = false;
    if (!cb.release())
      success = false;
    if (!h.release())
      success = false;
    if (!i1.release())
      success = false;
    if (!xi1.release())
      success = false;
    if (!i2.release())
      success = false;
    if (!xi2.release())
      success = false;
    if (!staging.release())
      success = false;
    device_ready = false;
    max_grid_x = 0;
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

__global__ void rsbb1e_cuda_kernel(double *sb, const double *cb, const double *h, const std::int64_t *i1, const double *xi1,
                                   const std::int64_t *i2, const double *xi2, std::size_t nrow, std::size_t nkastr,
                                   std::size_t nkaeff, std::size_t d2, std::size_t contribution_count) {
  const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (index >= contribution_count) {
    return;
  }

  const std::size_t row = index % nrow;
  const std::size_t traversal = index / nrow;
  const std::size_t q = traversal % nkaeff;
  const std::size_t i = traversal / nkaeff;
  const std::size_t i2_offset = q + i * nkastr;
  const std::int64_t dst = i2[i2_offset];
  if (dst == 0) {
    return;
  }

  double sum = 0.0;
  for (std::size_t j = 0; j < d2; ++j) {
    const std::size_t map_offset = q + j * nkastr;
    const std::int64_t src = i1[map_offset];
    if (src != 0) {
      sum += xi1[map_offset] * cb[row + static_cast<std::size_t>(src - 1) * nrow] * h[j + i * d2];
    }
  }
  atomic_add(sb + row + static_cast<std::size_t>(dst - 1) * nrow, xi2[i2_offset] * sum);
}

} // namespace

extern "C" int64_t lucia_rsbb1e_cuda_begin(double *sb, const double *cb, int64_t nrow) {
  if (sb == nullptr || cb == nullptr || workspace.session_active) {
    return 0;
  }
  std::size_t nrow_size = 0;
  if (!lucia_cuda::positive_size(nrow, &nrow_size) || !workspace.begin_session(sb, cb, nrow_size)) {
    return 0;
  }
  return 1;
}

extern "C" int64_t lucia_rsbb1e_cuda_end() {
  return workspace.end_session() ? 1 : -1;
}

extern "C" int64_t lucia_rsbb1e_cuda_route(double *sb, const double *cb, const double *h, const int64_t *i1, const double *xi1,
                                           const int64_t *i2, const double *xi2, int64_t nrow, int64_t ncb, int64_t nsb,
                                           int64_t nkastr, int64_t nkaeff, int64_t d1, int64_t d2, int64_t maxk) {
  if (sb == nullptr || cb == nullptr || h == nullptr || i1 == nullptr || xi1 == nullptr || i2 == nullptr || xi2 == nullptr) {
    return route_fallback(workspace, sb, cb);
  }
  if (d1 > 32 || d2 > 32 || nkaeff > nkastr) {
    return route_fallback(workspace, sb, cb);
  }

  std::size_t nrow_size = 0;
  std::size_t ncb_size = 0;
  std::size_t nsb_size = 0;
  std::size_t nkastr_size = 0;
  std::size_t nkaeff_size = 0;
  std::size_t d1_size = 0;
  std::size_t d2_size = 0;
  std::size_t maxk_size = 0;
  if (!lucia_cuda::positive_size(nrow, &nrow_size) || !lucia_cuda::positive_size(ncb, &ncb_size)
      || !lucia_cuda::positive_size(nsb, &nsb_size) || !lucia_cuda::positive_size(nkastr, &nkastr_size)
      || !lucia_cuda::positive_size(nkaeff, &nkaeff_size) || !lucia_cuda::positive_size(d1, &d1_size)
      || !lucia_cuda::positive_size(d2, &d2_size) || !lucia_cuda::positive_size(maxk, &maxk_size)) {
    return route_fallback(workspace, sb, cb);
  }

  std::size_t cb_count = 0;
  std::size_t sb_count = 0;
  std::size_t h_count = 0;
  std::size_t i1_count = 0;
  std::size_t xi1_count = 0;
  std::size_t i2_count = 0;
  std::size_t xi2_count = 0;
  if (!lucia_cuda::checked_mul(nrow_size, ncb_size, &cb_count) || !lucia_cuda::checked_mul(nrow_size, nsb_size, &sb_count)
      || !lucia_cuda::checked_mul(d2_size, d1_size, &h_count) || !lucia_cuda::checked_mul(nkastr_size, d2_size, &i1_count)
      || !lucia_cuda::checked_mul(nkastr_size, d2_size, &xi1_count) || !lucia_cuda::checked_mul(nkastr_size, d1_size, &i2_count)
      || !lucia_cuda::checked_mul(nkastr_size, d1_size, &xi2_count)) {
    return route_fallback(workspace, sb, cb);
  }

  const bool shared_match = lucia_sigma_cuda_blocks::sigma_blocks_match(sb, cb, sb_count, cb_count);
  const bool local_match = workspace.session_matches(sb, cb, nrow_size);
  if (workspace.session_active && !local_match && !shared_match) {
    return route_fallback(workspace, sb, cb);
  }
  if (local_match && workspace.session_resident
      && (workspace.session_sb_count != sb_count || workspace.session_cb_count != cb_count) && !workspace.flush_session()) {
    return -1;
  }

  const std::size_t nkaeff_scan = nkaeff_size;
  for (std::size_t j = 0; j < d2_size; ++j) {
    const std::size_t map_base = j * nkastr_size;
    for (std::size_t q = 0; q < nkaeff_scan; ++q) {
      const std::int64_t src = i1[map_base + q];
      if (src < 0 || src > ncb) {
        return route_fallback(workspace, sb, cb);
      }
    }
  }
  for (std::size_t i = 0; i < d1_size; ++i) {
    const std::size_t map_base = i * nkastr_size;
    for (std::size_t q = 0; q < nkaeff_scan; ++q) {
      const std::int64_t dst = i2[map_base + q];
      if (dst < 0 || dst > nsb) {
        return route_fallback(workspace, sb, cb);
      }
    }
  }

  std::size_t traversal_work = 0;
  if (!lucia_cuda::checked_mul(nrow_size, nkaeff_size, &traversal_work)
      || !lucia_cuda::checked_mul(traversal_work, d1_size, &traversal_work)
      || !lucia_cuda::checked_mul(traversal_work, d2_size, &traversal_work)) {
    return route_fallback(workspace, sb, cb);
  }

  std::size_t transfer_elements = 0;
  if (!lucia_cuda::checked_mul(sb_count, 2, &transfer_elements) || !checked_add(transfer_elements, cb_count, &transfer_elements)
      || !checked_add(transfer_elements, h_count, &transfer_elements)
      || !checked_add(transfer_elements, i1_count, &transfer_elements)
      || !checked_add(transfer_elements, xi1_count, &transfer_elements)
      || !checked_add(transfer_elements, i2_count, &transfer_elements)
      || !checked_add(transfer_elements, xi2_count, &transfer_elements)) {
    return route_fallback(workspace, sb, cb);
  }
  std::size_t transfer_bytes = 0;
  if (!lucia_cuda::checked_mul(transfer_elements, sizeof(double), &transfer_bytes)) {
    return route_fallback(workspace, sb, cb);
  }
  std::size_t minimum_transfer_work = 0;
  if (!lucia_cuda::checked_mul(transfer_bytes, minimum_work_per_transfer_byte, &minimum_transfer_work)) {
    return route_fallback(workspace, sb, cb);
  }
  if (traversal_work < minimum_traversal_work || traversal_work < minimum_transfer_work) {
    return route_fallback(workspace, sb, cb);
  }

  std::size_t contribution_count = 0;
  if (!lucia_cuda::checked_mul(nrow_size, nkaeff_size, &contribution_count)
      || !lucia_cuda::checked_mul(contribution_count, d1_size, &contribution_count)) {
    return route_fallback(workspace, sb, cb);
  }
  std::size_t grid_count = contribution_count / threads_per_block;
  if (contribution_count % threads_per_block != 0 && !checked_add(grid_count, 1, &grid_count)) {
    return route_fallback(workspace, sb, cb);
  }
  if (!workspace.supports_grid(grid_count)) {
    return route_fallback(workspace, sb, cb);
  }

  std::size_t sb_bytes = 0;
  if (!lucia_cuda::checked_mul(sb_count, sizeof(double), &sb_bytes)) {
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
  if (!workspace.h.copy_from(h, h_count) || !workspace.i1.copy_from(i1, i1_count) || !workspace.xi1.copy_from(xi1, xi1_count)
      || !workspace.i2.copy_from(i2, i2_count) || !workspace.xi2.copy_from(xi2, xi2_count)
      || (!shared_match && !local_match && !workspace.staging.reserve(sb_count))) {
    return route_fallback(workspace, sb, cb);
  }

  const dim3 block(static_cast<unsigned int>(threads_per_block), 1, 1);
  const dim3 grid(static_cast<unsigned int>(grid_count), 1, 1);
  rsbb1e_cuda_kernel<<<grid, block>>>(device_sb, device_cb, workspace.h.data, workspace.i1.data, workspace.xi1.data,
                                      workspace.i2.data, workspace.xi2.data, nrow_size, nkastr_size, nkaeff_size, d2_size,
                                      contribution_count);

  bool success = cudaGetLastError() == cudaSuccess;
  if (!shared_match && !local_match && cudaDeviceSynchronize() != cudaSuccess) {
    success = false;
  }
  if (!success) {
    workspace.invalidate_residency();
    return -1;
  }
  if (shared_match) {
    return 1;
  }
  if (local_match) {
    workspace.session_resident = true;
    workspace.session_sb_count = sb_count;
    workspace.session_cb_count = cb_count;
    return 1;
  }
  if (cudaMemcpy(workspace.staging.data, device_sb, sb_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) {
    workspace.invalidate_residency();
    return -1;
  }

  std::memcpy(sb, workspace.staging.data, sb_bytes);
  return 1;
}

extern "C" void lucia_rsbb1e_cuda_release() {
  workspace.release();
}
