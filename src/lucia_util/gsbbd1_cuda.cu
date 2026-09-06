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

struct Workspace {
  lucia_cuda::DeviceBuffer<double> rho1;
  lucia_cuda::DeviceBuffer<double> srho1;
  lucia_cuda::DeviceBuffer<double> sb;
  lucia_cuda::DeviceBuffer<double> cb;
  lucia_cuda::DeviceBuffer<double> xi1s;
  lucia_cuda::DeviceBuffer<double> xi2s;
  lucia_cuda::DeviceBuffer<std::int64_t> i1;
  lucia_cuda::DeviceBuffer<std::int64_t> i2;
  lucia_cuda::PinnedBuffer<double> rho1_staging;
  lucia_cuda::PinnedBuffer<double> srho1_staging;
  bool session_active = false;
  bool session_resident = false;
  double *session_rho1_host = nullptr;
  double *session_srho1_host = nullptr;
  const double *session_sb_host = nullptr;
  const double *session_cb_host = nullptr;
  std::size_t session_density_count = 0;
  std::size_t session_sb_count = 0;
  std::size_t session_cb_count = 0;
  bool session_update_spin = false;
  bool map_session_active = false;
  bool map_session_resident = false;
  const std::int64_t *map_i1_host = nullptr;
  const double *map_xi1s_host = nullptr;
  const std::int64_t *map_i2_host = nullptr;
  const double *map_xi2s_host = nullptr;
  std::size_t map_i1_count = 0;
  std::size_t map_i2_count = 0;

  bool session_matches(double *rho1_host, double *srho1_host, const double *sb_host, const double *cb_host,
                       std::size_t density_count, std::size_t sb_count, std::size_t cb_count, bool update_spin) const noexcept {
    return session_active && session_rho1_host == rho1_host && session_srho1_host == srho1_host && session_sb_host == sb_host
           && session_cb_host == cb_host && session_density_count == density_count && session_sb_count == sb_count
           && session_cb_count == cb_count && session_update_spin == update_spin;
  }

  bool begin_session(double *rho1_host, double *srho1_host, const double *sb_host, const double *cb_host, std::size_t density_count,
                     std::size_t sb_count, std::size_t cb_count, bool update_spin) noexcept {
    if (session_active) {
      return false;
    }
    session_active = true;
    session_resident = false;
    session_rho1_host = rho1_host;
    session_srho1_host = srho1_host;
    session_sb_host = sb_host;
    session_cb_host = cb_host;
    session_density_count = density_count;
    session_sb_count = sb_count;
    session_cb_count = cb_count;
    session_update_spin = update_spin;
    return true;
  }

  bool begin_map_session(const std::int64_t *i1_host, const double *xi1s_host, const std::int64_t *i2_host, const double *xi2s_host,
                         std::size_t i1_count, std::size_t i2_count) noexcept {
    if (map_session_active) {
      return false;
    }
    map_session_active = true;
    map_session_resident = false;
    map_i1_host = i1_host;
    map_xi1s_host = xi1s_host;
    map_i2_host = i2_host;
    map_xi2s_host = xi2s_host;
    map_i1_count = i1_count;
    map_i2_count = i2_count;
    return true;
  }

  bool map_session_matches(const std::int64_t *i1_host, const double *xi1s_host, const std::int64_t *i2_host,
                           const double *xi2s_host, std::size_t i1_count, std::size_t i2_count) const noexcept {
    return map_session_active && map_i1_host == i1_host && map_xi1s_host == xi1s_host && map_i2_host == i2_host
           && map_xi2s_host == xi2s_host && map_i1_count == i1_count && map_i2_count == i2_count;
  }

  void invalidate_map_residency() noexcept {
    map_session_resident = false;
  }

  void clear_map_session() noexcept {
    map_session_active = false;
    map_session_resident = false;
    map_i1_host = nullptr;
    map_xi1s_host = nullptr;
    map_i2_host = nullptr;
    map_xi2s_host = nullptr;
    map_i1_count = 0;
    map_i2_count = 0;
  }

  bool end_map_session() noexcept {
    clear_map_session();
    return true;
  }

  bool flush_session() noexcept {
    if (!session_active || !session_resident) {
      return true;
    }
    invalidate_residency();
    std::size_t density_bytes = 0;
    if (!lucia_cuda::checked_mul(session_density_count, sizeof(double), &density_bytes)
        || !rho1_staging.reserve(session_density_count) || (session_update_spin && !srho1_staging.reserve(session_density_count))) {
      return false;
    }
    const bool rho1_success = cudaMemcpy(rho1_staging.data, rho1.data, density_bytes, cudaMemcpyDeviceToHost) == cudaSuccess;
    const bool srho1_success
        = !session_update_spin || cudaMemcpy(srho1_staging.data, srho1.data, density_bytes, cudaMemcpyDeviceToHost) == cudaSuccess;
    if (!rho1_success || !srho1_success) {
      return false;
    }
    std::memcpy(session_rho1_host, rho1_staging.data, density_bytes);
    if (session_update_spin) {
      std::memcpy(session_srho1_host, srho1_staging.data, density_bytes);
    }
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
    session_rho1_host = nullptr;
    session_srho1_host = nullptr;
    session_sb_host = nullptr;
    session_cb_host = nullptr;
    session_density_count = 0;
    session_sb_count = 0;
    session_cb_count = 0;
    session_update_spin = false;
  }

  bool end_session() noexcept {
    const bool success = flush_session();
    clear_session();
    return success;
  }

  bool release() noexcept {
    bool success = end_session();
    clear_map_session();
    if (!rho1.release())
      success = false;
    if (!srho1.release())
      success = false;
    if (!sb.release())
      success = false;
    if (!cb.release())
      success = false;
    if (!xi1s.release())
      success = false;
    if (!xi2s.release())
      success = false;
    if (!i1.release())
      success = false;
    if (!i2.release())
      success = false;
    if (!rho1_staging.release())
      success = false;
    if (!srho1_staging.release())
      success = false;
    return success;
  }
};

__global__ void gsbbd1_cuda_kernel(double *rho1, double *srho1, const double *sb, const double *cb, const std::int64_t *i1,
                                   const double *xi1s, const std::int64_t *i2, const double *xi2s, std::size_t nacob,
                                   std::size_t nrow, std::size_t ibot0, std::size_t nibtc, std::size_t nkastr, std::size_t kbot0,
                                   std::size_t lkabtc, std::size_t d1, std::size_t off1, std::size_t off2, bool update_spin,
                                   double xab) {
  __shared__ double partial[block_size];
  const std::size_t pair = blockIdx.x;
  const std::size_t i = pair % d1;
  const std::size_t j = pair / d1;
  const std::size_t reduction = nibtc * lkabtc;
  double sum = 0.0;

  for (std::size_t index = threadIdx.x; index < reduction; index += blockDim.x) {
    const std::size_t k = index / nibtc;
    const std::size_t row = index - k * nibtc;
    const std::size_t c_map_offset = kbot0 + k + j * nkastr;
    const std::size_t s_map_offset = kbot0 + k + i * nkastr;
    const std::int64_t c_source = i1[c_map_offset];
    const std::int64_t s_source = i2[s_map_offset];
    if (c_source != 0 && s_source != 0) {
      const std::size_t c_offset = ibot0 + row + static_cast<std::size_t>(c_source - 1) * nrow;
      const std::size_t s_offset = ibot0 + row + static_cast<std::size_t>(s_source - 1) * nrow;
      sum += xi2s[s_map_offset] * sb[s_offset] * xi1s[c_map_offset] * cb[c_offset];
    }
  }
  partial[threadIdx.x] = sum;
  __syncthreads();

  for (std::size_t stride = block_size / 2; stride > 0; stride /= 2) {
    if (threadIdx.x < stride) {
      partial[threadIdx.x] += partial[threadIdx.x + stride];
    }
    __syncthreads();
  }
  if (threadIdx.x == 0) {
    const std::size_t output = off1 + i + (off2 + j) * nacob;
    rho1[output] += partial[0];
    if (update_spin) {
      srho1[output] += xab * partial[0];
    }
  }
}
Workspace workspace;

} // namespace

extern "C" int64_t lucia_gsbbd1_cuda_begin(double *rho1, double *srho1, const double *sb, const double *cb, int64_t nacob,
                                           int64_t nrow, int64_t nsb, int64_t ncb, int64_t idosrho1) {
  if (rho1 == nullptr || srho1 == nullptr || sb == nullptr || cb == nullptr || workspace.session_active) {
    return 0;
  }
  std::size_t nacob_size = 0;
  std::size_t nrow_size = 0;
  std::size_t nsb_size = 0;
  std::size_t ncb_size = 0;
  std::size_t density_count = 0;
  std::size_t sb_count = 0;
  std::size_t cb_count = 0;
  if (!lucia_cuda::positive_size(nacob, &nacob_size) || !lucia_cuda::positive_size(nrow, &nrow_size)
      || !lucia_cuda::positive_size(nsb, &nsb_size) || !lucia_cuda::positive_size(ncb, &ncb_size) || idosrho1 < 0 || idosrho1 > 1
      || !lucia_cuda::checked_mul(nacob_size, nacob_size, &density_count)
      || !lucia_cuda::checked_mul(nrow_size, nsb_size, &sb_count) || !lucia_cuda::checked_mul(nrow_size, ncb_size, &cb_count)
      || !workspace.begin_session(rho1, srho1, sb, cb, density_count, sb_count, cb_count, idosrho1 == 1)) {
    return 0;
  }
  return 1;
}

extern "C" int64_t lucia_gsbbd1_cuda_maps_begin(const int64_t *i1, const double *xi1s, const int64_t *i2, const double *xi2s,
                                                int64_t nkastr, int64_t d1, int64_t d2) {
  if (i1 == nullptr || xi1s == nullptr || i2 == nullptr || xi2s == nullptr) {
    return 0;
  }
  std::size_t nkastr_size = 0;
  std::size_t d1_size = 0;
  std::size_t d2_size = 0;
  std::size_t i1_count = 0;
  std::size_t i2_count = 0;
  if (!lucia_cuda::positive_size(nkastr, &nkastr_size) || !lucia_cuda::positive_size(d1, &d1_size)
      || !lucia_cuda::positive_size(d2, &d2_size) || !lucia_cuda::checked_mul(nkastr_size, d2_size, &i1_count)
      || !lucia_cuda::checked_mul(nkastr_size, d1_size, &i2_count)
      || !workspace.begin_map_session(i1, xi1s, i2, xi2s, i1_count, i2_count)) {
    return 0;
  }
  return 1;
}

extern "C" int64_t lucia_gsbbd1_cuda_maps_end() {
  return workspace.end_map_session() ? 1 : -1;
}

extern "C" int64_t lucia_gsbbd1_cuda_end() {
  return workspace.end_session() ? 1 : -1;
}

extern "C" int64_t lucia_gsbbd1_cuda_route(double *rho1, double *srho1, const double *sb, const double *cb, const int64_t *i1,
                                           const double *xi1s, const int64_t *i2, const double *xi2s, int64_t nacob, int64_t nrow,
                                           int64_t nsb, int64_t ncb, int64_t ibot, int64_t nibtc, int64_t nkastr, int64_t kbot,
                                           int64_t lkabtc, int64_t d1, int64_t d2, int64_t off1, int64_t off2, int64_t idosrho1,
                                           double xab) {
  if (rho1 == nullptr || srho1 == nullptr || sb == nullptr || cb == nullptr || i1 == nullptr || xi1s == nullptr || i2 == nullptr
      || xi2s == nullptr || !std::isfinite(xab)) {
    return workspace.fallback_status();
  }

  std::size_t nacob_size = 0;
  std::size_t nrow_size = 0;
  std::size_t nsb_size = 0;
  std::size_t ncb_size = 0;
  std::size_t ibot_size = 0;
  std::size_t nibtc_size = 0;
  std::size_t nkastr_size = 0;
  std::size_t kbot_size = 0;
  std::size_t lkabtc_size = 0;
  std::size_t d1_size = 0;
  std::size_t d2_size = 0;
  if (!lucia_cuda::positive_size(nacob, &nacob_size) || !lucia_cuda::positive_size(nrow, &nrow_size)
      || !lucia_cuda::positive_size(nsb, &nsb_size) || !lucia_cuda::positive_size(ncb, &ncb_size)
      || !lucia_cuda::positive_size(ibot, &ibot_size) || !lucia_cuda::positive_size(nibtc, &nibtc_size)
      || !lucia_cuda::positive_size(nkastr, &nkastr_size) || !lucia_cuda::positive_size(kbot, &kbot_size)
      || !lucia_cuda::positive_size(lkabtc, &lkabtc_size) || !lucia_cuda::positive_size(d1, &d1_size)
      || !lucia_cuda::positive_size(d2, &d2_size) || off1 < 0 || off2 < 0 || idosrho1 < 0 || idosrho1 > 1) {
    return workspace.fallback_status();
  }

  const std::size_t ibot0 = ibot_size - 1;
  const std::size_t kbot0 = kbot_size - 1;
  const std::size_t off1_size = static_cast<std::size_t>(off1);
  const std::size_t off2_size = static_cast<std::size_t>(off2);
  if (ibot0 >= nrow_size || nibtc_size > nrow_size - ibot0 || kbot0 >= nkastr_size || lkabtc_size > nkastr_size - kbot0
      || off1_size >= nacob_size || d1_size > nacob_size - off1_size || off2_size >= nacob_size
      || d2_size > nacob_size - off2_size) {
    return workspace.fallback_status();
  }

  std::size_t density_count = 0;
  std::size_t sb_count = 0;
  std::size_t cb_count = 0;
  std::size_t i1_count = 0;
  std::size_t i2_count = 0;
  std::size_t pair_count = 0;
  if (!lucia_cuda::checked_mul(nacob_size, nacob_size, &density_count) || !lucia_cuda::checked_mul(nrow_size, nsb_size, &sb_count)
      || !lucia_cuda::checked_mul(nrow_size, ncb_size, &cb_count) || !lucia_cuda::checked_mul(nkastr_size, d2_size, &i1_count)
      || !lucia_cuda::checked_mul(nkastr_size, d1_size, &i2_count) || !lucia_cuda::checked_mul(d1_size, d2_size, &pair_count)
      || pair_count > static_cast<std::size_t>((std::numeric_limits<unsigned int>::max)())) {
    return workspace.fallback_status();
  }

  const bool map_match = workspace.map_session_matches(i1, xi1s, i2, xi2s, i1_count, i2_count);
  if (workspace.map_session_active && !map_match) {
    return workspace.fallback_status();
  }

  const bool update_spin = idosrho1 == 1;
  const bool session_match = workspace.session_matches(rho1, srho1, sb, cb, density_count, sb_count, cb_count, update_spin);
  if (workspace.session_active && !session_match && workspace.fallback_status() != 0) {
    return -1;
  }
  const bool shared_match = lucia_sigma_cuda_blocks::sigma_blocks_match(const_cast<double *>(sb), cb, sb_count, cb_count);

  for (std::size_t j = 0; j < d2_size; ++j) {
    for (std::size_t k = 0; k < lkabtc_size; ++k) {
      const std::size_t offset = kbot0 + k + j * nkastr_size;
      if (i1[offset] < 0 || i1[offset] > ncb || !std::isfinite(xi1s[offset])) {
        return workspace.fallback_status();
      }
    }
  }
  for (std::size_t i = 0; i < d1_size; ++i) {
    for (std::size_t k = 0; k < lkabtc_size; ++k) {
      const std::size_t offset = kbot0 + k + i * nkastr_size;
      if (i2[offset] < 0 || i2[offset] > nsb || !std::isfinite(xi2s[offset])) {
        return workspace.fallback_status();
      }
    }
  }

  const bool upload_blocks = session_match && !workspace.session_resident;
  const bool upload_maps = map_match && !workspace.map_session_resident;
  double *device_sb = workspace.sb.data;
  const double *device_cb = workspace.cb.data;
  if (shared_match
      && !lucia_sigma_cuda_blocks::sigma_blocks_acquire(const_cast<double *>(sb), cb, sb_count, cb_count, &device_sb, &device_cb)) {
    return workspace.fallback_status();
  }
  if ((!session_match || upload_blocks)
      && (!workspace.rho1.copy_from(rho1, density_count) || (update_spin && !workspace.srho1.copy_from(srho1, density_count)))) {
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
  if ((!map_match || upload_maps)
      && (!workspace.i1.copy_from(i1, i1_count) || !workspace.xi1s.copy_from(xi1s, i1_count)
          || !workspace.i2.copy_from(i2, i2_count) || !workspace.xi2s.copy_from(xi2s, i2_count))) {
    return workspace.fallback_status();
  }
  if (map_match && upload_maps) {
    workspace.map_session_resident = true;
  }
  if ((!session_match || upload_blocks)
      && (!workspace.rho1_staging.reserve(density_count) || (update_spin && !workspace.srho1_staging.reserve(density_count)))) {
    return workspace.fallback_status();
  }

  gsbbd1_cuda_kernel<<<static_cast<unsigned int>(pair_count), block_size>>>(
      workspace.rho1.data, workspace.srho1.data, device_sb, device_cb, workspace.i1.data, workspace.xi1s.data, workspace.i2.data,
      workspace.xi2s.data, nacob_size, nrow_size, ibot0, nibtc_size, nkastr_size, kbot0, lkabtc_size, d1_size, off1_size, off2_size,
      update_spin, xab);
  bool success = cudaGetLastError() == cudaSuccess;
  if (!session_match && cudaDeviceSynchronize() != cudaSuccess) {
    success = false;
  }

  if (!success) {
    workspace.invalidate_residency();
    workspace.invalidate_map_residency();
    return -1;
  }
  if (session_match) {
    workspace.session_resident = true;
    return 1;
  }

  std::size_t density_bytes = 0;
  if (!lucia_cuda::checked_mul(density_count, sizeof(double), &density_bytes)
      || cudaMemcpy(workspace.rho1_staging.data, workspace.rho1.data, density_bytes, cudaMemcpyDeviceToHost) != cudaSuccess
      || (update_spin
          && cudaMemcpy(workspace.srho1_staging.data, workspace.srho1.data, density_bytes, cudaMemcpyDeviceToHost)
                 != cudaSuccess)) {
    workspace.invalidate_residency();
    return -1;
  }
  std::memcpy(rho1, workspace.rho1_staging.data, density_bytes);
  if (update_spin) {
    std::memcpy(srho1, workspace.srho1_staging.data, density_bytes);
  }
  return 1;
}

extern "C" void lucia_gsbbd1_cuda_release() {
  workspace.release();
}
