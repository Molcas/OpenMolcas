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

bool checked_add(std::size_t left, std::size_t right, std::size_t *result) noexcept
{
  if (left > (std::numeric_limits<std::size_t>::max)() - right) {
    return false;
  }
  *result = left + right;
  return true;
}
bool checked_product3(std::size_t a, std::size_t b, std::size_t c, std::size_t *result) noexcept
{
  std::size_t partial = 0;
  return lucia_cuda::checked_mul(a, b, &partial) && lucia_cuda::checked_mul(partial, c, result);
}

bool checked_product4(std::size_t a, std::size_t b, std::size_t c, std::size_t d, std::size_t *result) noexcept
{
  std::size_t partial = 0;
  return checked_product3(a, b, c, &partial) && lucia_cuda::checked_mul(partial, d, result);
}

bool append_section(std::size_t bytes, std::size_t *offset, std::size_t *total) noexcept
{
  *offset = *total;
  return checked_add(*total, bytes, total);
}

struct RouteSlot {
  lucia_cuda::DeviceBuffer<unsigned char> device;
  lucia_cuda::PinnedBuffer<unsigned char> host;
  cudaEvent_t host_ready = nullptr;
  bool pending = false;

  bool prepare() noexcept
  {
    if (!pending) {
      return true;
    }
    if (cudaEventSynchronize(host_ready) != cudaSuccess) {
      return false;
    }
    pending = false;
    return true;
  }

  bool ensure_event() noexcept
  {
    return host_ready != nullptr ||
           cudaEventCreateWithFlags(&host_ready, cudaEventDisableTiming) == cudaSuccess;
  }

  bool upload(std::size_t bytes) noexcept
  {
    if (!device.reserve(bytes) || !ensure_event() ||
        cudaMemcpyAsync(device.data, host.data, bytes, cudaMemcpyHostToDevice, 0) != cudaSuccess ||
        cudaEventRecord(host_ready, 0) != cudaSuccess) {
      return false;
    }
    pending = true;
    return true;
  }

  bool release() noexcept
  {
    bool success = true;
    if (pending && cudaEventSynchronize(host_ready) != cudaSuccess) {
      success = false;
    }
    pending = false;
    if (host_ready != nullptr && cudaEventDestroy(host_ready) != cudaSuccess) {
      success = false;
    }
    host_ready = nullptr;
    if (!device.release()) success = false;
    if (!host.release()) success = false;
    return success;
  }
};

struct Workspace {
  lucia_cuda::DeviceBuffer<double> sb;
  lucia_cuda::DeviceBuffer<double> cb;
  RouteSlot route[2];
  unsigned int next_route = 0;
  lucia_cuda::PinnedBuffer<std::int32_t> inverse_cursor_host;
  lucia_cuda::PinnedBuffer<double> staging;
  bool device_ready = false;
  bool session_active = false;
  bool shared_session_active = false;
  bool session_resident = false;
  double *session_sb_host = nullptr;
  const double *session_cb_host = nullptr;
  std::size_t session_sb_count = 0;
  std::size_t session_cb_count = 0;
  std::size_t max_grid_x = 0;

  bool session_matches(double *sb_host, const double *cb_host,
                       std::size_t sb_count, std::size_t cb_count) const noexcept
  {
    return session_active && session_sb_host == sb_host && session_cb_host == cb_host &&
           session_sb_count == sb_count && session_cb_count == cb_count;
  }

  bool begin_session(double *sb_host, const double *cb_host,
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

  bool flush_session() noexcept
  {
    if (!session_active || !session_resident) {
      return true;
    }
    invalidate_residency();
    std::size_t sb_bytes = 0;
    if (!lucia_cuda::checked_mul(session_sb_count, sizeof(double), &sb_bytes) ||
        !staging.reserve(session_sb_count) ||
        cudaMemcpy(staging.data, sb.data, sb_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) {
      return false;
    }
    std::memcpy(session_sb_host, staging.data, sb_bytes);
    return true;
  }

  void invalidate_residency() noexcept
  {
    session_resident = false;
  }

  std::int64_t fallback_status() noexcept
  {
    return flush_session() ? 0 : -1;
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
    const bool success = flush_session();
    clear_session();
    return success;
  }

  bool supports_launch(std::size_t block_count) noexcept
  {
    if (!device_ready) {
      int device = 0;
      cudaDeviceProp properties{};
      if (cudaGetDevice(&device) != cudaSuccess ||
          cudaGetDeviceProperties(&properties, device) != cudaSuccess ||
          properties.maxGridSize[0] <= 0 || properties.maxThreadsPerBlock < 256) {
        return false;
      }
      max_grid_x = static_cast<std::size_t>(properties.maxGridSize[0]);
      device_ready = true;
    }
    return block_count > 0 && block_count <= max_grid_x &&
           block_count <= static_cast<std::size_t>((std::numeric_limits<unsigned int>::max)());
  }

  bool release() noexcept
  {
    bool success = end_session();
    shared_session_active = false;
    if (!sb.release()) success = false;
    if (!cb.release()) success = false;
    if (!route[0].release()) success = false;
    if (!route[1].release()) success = false;
    next_route = 0;
    if (!inverse_cursor_host.release()) success = false;
    if (!staging.release()) success = false;
    device_ready = false;
    max_grid_x = 0;
    return success;
  }
};

std::int64_t route_fallback(Workspace &workspace, double *sb, const double *cb) noexcept
{
  bool matched = false;
  if (!lucia_sigma_cuda_blocks::sigma_blocks_flush_if_match(sb, cb, &matched)) {
    return -1;
  }
  return matched ? 0 : workspace.fallback_status();
}

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

__global__ void rsbb2bn_cuda_kernel(
    double *sb, const double *cb, const double *xint,
    const std::int64_t *i3, const double *xi3s,
    const std::int32_t *inverse_offset, const std::int32_t *inverse_k,
    const std::int32_t *inverse_kb, const double *inverse_factor,
    const std::int32_t *l_offset, const std::int32_t *l_index,
    const std::int64_t *l_jb, const double *l_factor,
    const std::int32_t *j_offset, const std::int32_t *j_index,
    const std::int64_t *j_ja, const double *j_factor,
    std::size_t nib, std::size_t njb, std::size_t nkastr,
    std::size_t kabot0, std::size_t lkabtc, std::size_t ni,
    std::size_t nj, std::size_t nk, std::size_t nl, bool ordered)
{
  const std::size_t block = static_cast<std::size_t>(blockIdx.x);
  const std::size_t i = block / lkabtc;
  const std::size_t ka = block - i * lkabtc;
  const std::size_t scatter_map = kabot0 + ka + i * nkastr;
  const std::int64_t ia = i3[scatter_map];
  if (ia == 0) {
    return;
  }

  for (std::size_t ib = threadIdx.x; ib < nib; ib += blockDim.x) {
    double result = 0.0;
    for (std::int32_t p = inverse_offset[ib]; p < inverse_offset[ib + 1]; ++p) {
      const std::size_t k = static_cast<std::size_t>(inverse_k[p]);
      if (ordered && i > k) {
        continue;
      }
      const std::size_t kb = static_cast<std::size_t>(inverse_kb[p]);
      double kl_sum = 0.0;
      for (std::int32_t lp = l_offset[kb]; lp < l_offset[kb + 1]; ++lp) {
        const std::size_t l = static_cast<std::size_t>(l_index[lp]);
        double sum = 0.0;
        for (std::int32_t jp = j_offset[ka]; jp < j_offset[ka + 1]; ++jp) {
          const std::size_t j = static_cast<std::size_t>(j_index[jp]);
          double integral = xint[j + i * nj + (l * nk + k) * ni * nj];
          if (ordered && i == k) {
            if (j > l && j < nl) {
              integral = 0.0;
            } else if (j == l) {
              integral *= 0.5;
            }
          }
          const std::size_t c_offset = static_cast<std::size_t>(l_jb[lp] - 1) +
                                       static_cast<std::size_t>(j_ja[jp] - 1) * njb;
          sum += j_factor[jp] * cb[c_offset] * integral;
        }
        kl_sum += l_factor[lp] * sum;
      }
      result += inverse_factor[p] * kl_sum;
    }
    if (result != 0.0) {
      atomic_add(sb + ib + static_cast<std::size_t>(ia - 1) * nib,
                 xi3s[scatter_map] * result);
    }
  }
}

Workspace workspace;

} // namespace

extern "C" int64_t lucia_rsbb2bn_cuda_begin(
    double *sb, const double *cb, int64_t nia, int64_t nib, int64_t nja, int64_t njb)
{
  if (sb == nullptr || cb == nullptr) {
    return 0;
  }
  std::size_t nia_size = 0;
  std::size_t nib_size = 0;
  std::size_t nja_size = 0;
  std::size_t njb_size = 0;
  std::size_t sb_count = 0;
  std::size_t cb_count = 0;
  if (!lucia_cuda::positive_size(nia, &nia_size) ||
      !lucia_cuda::positive_size(nib, &nib_size) ||
      !lucia_cuda::positive_size(nja, &nja_size) ||
      !lucia_cuda::positive_size(njb, &njb_size) ||
      !lucia_cuda::checked_mul(nia_size, nib_size, &sb_count) ||
      !lucia_cuda::checked_mul(nja_size, njb_size, &cb_count)) {
    return 0;
  }
  if (lucia_sigma_cuda_blocks::sigma_blocks_match(sb, cb, sb_count, cb_count)) {
    workspace.shared_session_active = true;
    return 1;
  }
  if (workspace.session_active || !workspace.begin_session(sb, cb, sb_count, cb_count)) {
    return 0;
  }
  workspace.shared_session_active = false;
  return 1;
}

extern "C" int64_t lucia_rsbb2bn_cuda_flush()
{
  if (workspace.shared_session_active) {
    bool active = false;
    if (!lucia_sigma_cuda_blocks::sigma_blocks_flush_active(&active)) {
      return -1;
    }
    return active ? 1 : -1;
  }
  return workspace.flush_session() ? 1 : -1;
}

extern "C" int64_t lucia_rsbb2bn_cuda_end()
{
  if (workspace.shared_session_active) {
    workspace.shared_session_active = false;
    return 1;
  }
  return workspace.end_session() ? 1 : -1;
}

extern "C" int64_t lucia_rsbb2bn_cuda_route(
    double *sb, const double *cb, const double *xint, const int64_t *i1,
    const double *xi1s, const int64_t *i3, const double *xi3s, const int64_t *i4,
    const double *xi4s, const int64_t *i2, const double *xi2s, int64_t nia,
    int64_t nib, int64_t nja, int64_t njb, int64_t nkastr, int64_t kabot,
    int64_t lkabtc, int64_t nkbstr, int64_t ni, int64_t nj, int64_t nk,
    int64_t nl, int64_t ikord)
{
  if (sb == nullptr || cb == nullptr || xint == nullptr || i1 == nullptr || xi1s == nullptr ||
      i3 == nullptr || xi3s == nullptr || i4 == nullptr || xi4s == nullptr || i2 == nullptr ||
      xi2s == nullptr || (ikord != 0 && ikord != 1)) {
    return route_fallback(workspace, sb, cb);
  }

  std::size_t nia_size = 0, nib_size = 0, nja_size = 0, njb_size = 0;
  std::size_t nkastr_size = 0, kabot_size = 0, lkabtc_size = 0, nkbstr_size = 0;
  std::size_t ni_size = 0, nj_size = 0, nk_size = 0, nl_size = 0;
  if (!lucia_cuda::positive_size(nia, &nia_size) || !lucia_cuda::positive_size(nib, &nib_size) ||
      !lucia_cuda::positive_size(nja, &nja_size) || !lucia_cuda::positive_size(njb, &njb_size) ||
      !lucia_cuda::positive_size(nkastr, &nkastr_size) || !lucia_cuda::positive_size(kabot, &kabot_size) ||
      !lucia_cuda::positive_size(lkabtc, &lkabtc_size) || !lucia_cuda::positive_size(nkbstr, &nkbstr_size) ||
      !lucia_cuda::positive_size(ni, &ni_size) || !lucia_cuda::positive_size(nj, &nj_size) ||
      !lucia_cuda::positive_size(nk, &nk_size) || !lucia_cuda::positive_size(nl, &nl_size)) {
    return route_fallback(workspace, sb, cb);
  }

  const std::size_t kabot0 = kabot_size - 1;
  if (kabot0 >= nkastr_size || lkabtc_size > nkastr_size - kabot0 ||
      (ikord == 1 && (nk_size > ni_size || nl_size > nj_size))) {
    return route_fallback(workspace, sb, cb);
  }

  std::size_t sb_count = 0, cb_count = 0, xint_count = 0;
  std::size_t i1_count = 0, i3_count = 0, i4_count = 0, i2_count = 0;
  std::size_t operation_count = 0, work_count = 0, block_count = 0;
  if (!lucia_cuda::checked_mul(nia_size, nib_size, &sb_count) ||
      !lucia_cuda::checked_mul(nja_size, njb_size, &cb_count) ||
      !checked_product4(ni_size, nj_size, nk_size, nl_size, &xint_count) ||
      !lucia_cuda::checked_mul(nkastr_size, nj_size, &i1_count) ||
      !lucia_cuda::checked_mul(nkastr_size, ni_size, &i3_count) ||
      !lucia_cuda::checked_mul(nkbstr_size, nk_size, &i4_count) ||
      !lucia_cuda::checked_mul(nkbstr_size, nl_size, &i2_count) ||
      !checked_product3(nkbstr_size, nk_size, nl_size, &operation_count) ||
      !checked_product4(lkabtc_size, ni_size, nj_size, operation_count, &work_count) ||
      !lucia_cuda::checked_mul(lkabtc_size, ni_size, &block_count)) {
    return route_fallback(workspace, sb, cb);
  }

  const bool shared_match = lucia_sigma_cuda_blocks::sigma_blocks_match(sb, cb, sb_count, cb_count);
  const bool local_match = workspace.session_matches(sb, cb, sb_count, cb_count);
  if (workspace.session_active && !local_match && !shared_match &&
      route_fallback(workspace, sb, cb) != 0) {
    return -1;
  }

  constexpr std::size_t minimum_work = 300000000;
  if (work_count < minimum_work || !workspace.supports_launch(block_count)) {
    return route_fallback(workspace, sb, cb);
  }

  for (std::size_t j = 0; j < nj_size; ++j) {
    for (std::size_t ka = 0; ka < lkabtc_size; ++ka) {
      const std::size_t offset = kabot0 + ka + j * nkastr_size;
      if (i1[offset] < 0 || i1[offset] > nja || !std::isfinite(xi1s[offset])) {
        return route_fallback(workspace, sb, cb);
      }
    }
  }
  for (std::size_t i = 0; i < ni_size; ++i) {
    for (std::size_t ka = 0; ka < lkabtc_size; ++ka) {
      const std::size_t offset = kabot0 + ka + i * nkastr_size;
      if (i3[offset] < 0 || i3[offset] > nia || !std::isfinite(xi3s[offset])) {
        return route_fallback(workspace, sb, cb);
      }
    }
  }
  for (std::size_t offset = 0; offset < i4_count; ++offset) {
    if (i4[offset] < 0 || i4[offset] > nib || !std::isfinite(xi4s[offset])) {
      return route_fallback(workspace, sb, cb);
    }
  }
  for (std::size_t offset = 0; offset < i2_count; ++offset) {
    if (i2[offset] < 0 || i2[offset] > njb || !std::isfinite(xi2s[offset])) {
      return route_fallback(workspace, sb, cb);
    }
  }

  std::size_t inverse_offset_count = 0, l_offset_count = 0, j_offset_count = 0;
  std::size_t j_capacity = 0;
  const std::size_t max_compact =
      static_cast<std::size_t>((std::numeric_limits<std::int32_t>::max)());
  if (i4_count > max_compact || i2_count > max_compact ||
      nk_size > max_compact || nkbstr_size > max_compact ||
      nl_size > max_compact || nj_size > max_compact ||
      !lucia_cuda::checked_mul(lkabtc_size, nj_size, &j_capacity) ||
      j_capacity > max_compact || !checked_add(nib_size, 1, &inverse_offset_count) ||
      !checked_add(nkbstr_size, 1, &l_offset_count) ||
      !checked_add(lkabtc_size, 1, &j_offset_count)) {
    return route_fallback(workspace, sb, cb);
  }

  std::size_t inverse_count = 0, l_count = 0, j_count = 0;
  for (std::size_t k = 0; k < nk_size; ++k) {
    for (std::size_t kb = 0; kb < nkbstr_size; ++kb) {
      const std::size_t map = kb + k * nkbstr_size;
      if (i4[map] != 0) {
        ++inverse_count;
      }
    }
  }
  for (std::size_t kb = 0; kb < nkbstr_size; ++kb) {
    for (std::size_t l = 0; l < nl_size; ++l) {
      if (i2[kb + l * nkbstr_size] != 0) {
        ++l_count;
      }
    }
  }
  for (std::size_t ka = 0; ka < lkabtc_size; ++ka) {
    for (std::size_t j = 0; j < nj_size; ++j) {
      if (i1[kabot0 + ka + j * nkastr_size] != 0) {
        ++j_count;
      }
    }
  }
  if (inverse_count > max_compact || l_count > max_compact || j_count > max_compact) {
    return route_fallback(workspace, sb, cb);
  }

  std::size_t xint_bytes = 0, i3_bytes = 0, xi3s_bytes = 0;
  std::size_t inverse_factor_bytes = 0, l_jb_bytes = 0, l_factor_bytes = 0;
  std::size_t j_ja_bytes = 0, j_factor_bytes = 0;
  std::size_t inverse_offset_bytes = 0, inverse_k_bytes = 0, inverse_kb_bytes = 0;
  std::size_t l_offset_bytes = 0, l_index_bytes = 0;
  std::size_t j_offset_bytes = 0, j_index_bytes = 0;
  if (!lucia_cuda::checked_mul(xint_count, sizeof(double), &xint_bytes) ||
      !lucia_cuda::checked_mul(i3_count, sizeof(std::int64_t), &i3_bytes) ||
      !lucia_cuda::checked_mul(i3_count, sizeof(double), &xi3s_bytes) ||
      !lucia_cuda::checked_mul(inverse_count, sizeof(double), &inverse_factor_bytes) ||
      !lucia_cuda::checked_mul(l_count, sizeof(std::int64_t), &l_jb_bytes) ||
      !lucia_cuda::checked_mul(l_count, sizeof(double), &l_factor_bytes) ||
      !lucia_cuda::checked_mul(j_count, sizeof(std::int64_t), &j_ja_bytes) ||
      !lucia_cuda::checked_mul(j_count, sizeof(double), &j_factor_bytes) ||
      !lucia_cuda::checked_mul(inverse_offset_count, sizeof(std::int32_t), &inverse_offset_bytes) ||
      !lucia_cuda::checked_mul(inverse_count, sizeof(std::int32_t), &inverse_k_bytes) ||
      !lucia_cuda::checked_mul(inverse_count, sizeof(std::int32_t), &inverse_kb_bytes) ||
      !lucia_cuda::checked_mul(l_offset_count, sizeof(std::int32_t), &l_offset_bytes) ||
      !lucia_cuda::checked_mul(l_count, sizeof(std::int32_t), &l_index_bytes) ||
      !lucia_cuda::checked_mul(j_offset_count, sizeof(std::int32_t), &j_offset_bytes) ||
      !lucia_cuda::checked_mul(j_count, sizeof(std::int32_t), &j_index_bytes)) {
    return route_fallback(workspace, sb, cb);
  }

  std::size_t inverse_factor_offset = 0, l_jb_offset = 0, l_factor_offset = 0;
  std::size_t j_ja_offset = 0, j_factor_offset = 0;
  std::size_t inverse_offset_offset = 0, inverse_k_offset = 0, inverse_kb_offset = 0;
  std::size_t l_offset_offset = 0, l_index_offset = 0;
  std::size_t j_offset_offset = 0, j_index_offset = 0;
  std::size_t xint_offset = 0, i3_offset = 0, xi3s_offset = 0;
  std::size_t payload_bytes = 0;
  if (!append_section(xint_bytes, &xint_offset, &payload_bytes) ||
      !append_section(i3_bytes, &i3_offset, &payload_bytes) ||
      !append_section(xi3s_bytes, &xi3s_offset, &payload_bytes) ||
      !append_section(inverse_factor_bytes, &inverse_factor_offset, &payload_bytes) ||
      !append_section(l_jb_bytes, &l_jb_offset, &payload_bytes) ||
      !append_section(l_factor_bytes, &l_factor_offset, &payload_bytes) ||
      !append_section(j_ja_bytes, &j_ja_offset, &payload_bytes) ||
      !append_section(j_factor_bytes, &j_factor_offset, &payload_bytes) ||
      !append_section(inverse_offset_bytes, &inverse_offset_offset, &payload_bytes) ||
      !append_section(inverse_k_bytes, &inverse_k_offset, &payload_bytes) ||
      !append_section(inverse_kb_bytes, &inverse_kb_offset, &payload_bytes) ||
      !append_section(l_offset_bytes, &l_offset_offset, &payload_bytes) ||
      !append_section(l_index_bytes, &l_index_offset, &payload_bytes) ||
      !append_section(j_offset_bytes, &j_offset_offset, &payload_bytes) ||
      !append_section(j_index_bytes, &j_index_offset, &payload_bytes)) {
    return route_fallback(workspace, sb, cb);
  }

  RouteSlot &slot = workspace.route[workspace.next_route];
  if (!slot.prepare() || !slot.host.reserve(payload_bytes) ||
      !workspace.inverse_cursor_host.reserve(nib_size)) {
    return route_fallback(workspace, sb, cb);
  }

  unsigned char *payload_host = slot.host.data;
  std::memcpy(payload_host + xint_offset, xint, xint_bytes);
  std::memcpy(payload_host + i3_offset, i3, i3_bytes);
  std::memcpy(payload_host + xi3s_offset, xi3s, xi3s_bytes);
  double *inverse_factor_host = reinterpret_cast<double *>(payload_host + inverse_factor_offset);
  std::int64_t *l_jb_host = reinterpret_cast<std::int64_t *>(payload_host + l_jb_offset);
  double *l_factor_host = reinterpret_cast<double *>(payload_host + l_factor_offset);
  std::int64_t *j_ja_host = reinterpret_cast<std::int64_t *>(payload_host + j_ja_offset);
  double *j_factor_host = reinterpret_cast<double *>(payload_host + j_factor_offset);
  std::int32_t *inverse_offset_host =
      reinterpret_cast<std::int32_t *>(payload_host + inverse_offset_offset);
  std::int32_t *inverse_k_host = reinterpret_cast<std::int32_t *>(payload_host + inverse_k_offset);
  std::int32_t *inverse_kb_host = reinterpret_cast<std::int32_t *>(payload_host + inverse_kb_offset);
  std::int32_t *l_offset_host = reinterpret_cast<std::int32_t *>(payload_host + l_offset_offset);
  std::int32_t *l_index_host = reinterpret_cast<std::int32_t *>(payload_host + l_index_offset);
  std::int32_t *j_offset_host = reinterpret_cast<std::int32_t *>(payload_host + j_offset_offset);
  std::int32_t *j_index_host = reinterpret_cast<std::int32_t *>(payload_host + j_index_offset);

  for (std::size_t ib = 0; ib < inverse_offset_count; ++ib) {
    inverse_offset_host[ib] = 0;
  }
  for (std::size_t k = 0; k < nk_size; ++k) {
    for (std::size_t kb = 0; kb < nkbstr_size; ++kb) {
      const std::size_t map = kb + k * nkbstr_size;
      if (i4[map] != 0) {
        ++inverse_offset_host[static_cast<std::size_t>(i4[map])];
      }
    }
  }
  for (std::size_t ib = 0; ib < nib_size; ++ib) {
    inverse_offset_host[ib + 1] += inverse_offset_host[ib];
    workspace.inverse_cursor_host.data[ib] = inverse_offset_host[ib];
  }
  std::size_t inverse_filled = 0;
  for (std::size_t k = 0; k < nk_size; ++k) {
    for (std::size_t kb = 0; kb < nkbstr_size; ++kb) {
      const std::size_t map = kb + k * nkbstr_size;
      if (i4[map] != 0) {
        const std::size_t ib = static_cast<std::size_t>(i4[map] - 1);
        const std::int32_t position = workspace.inverse_cursor_host.data[ib]++;
        inverse_k_host[position] = static_cast<std::int32_t>(k);
        inverse_kb_host[position] = static_cast<std::int32_t>(kb);
        inverse_factor_host[position] = xi4s[map];
        ++inverse_filled;
      }
    }
  }
  if (inverse_filled != inverse_count) {
    return route_fallback(workspace, sb, cb);
  }

  std::size_t l_filled = 0;
  for (std::size_t kb = 0; kb < nkbstr_size; ++kb) {
    l_offset_host[kb] = static_cast<std::int32_t>(l_filled);
    for (std::size_t l = 0; l < nl_size; ++l) {
      const std::size_t map = kb + l * nkbstr_size;
      if (i2[map] != 0) {
        l_index_host[l_filled] = static_cast<std::int32_t>(l);
        l_jb_host[l_filled] = i2[map];
        l_factor_host[l_filled] = xi2s[map];
        ++l_filled;
      }
    }
  }
  l_offset_host[nkbstr_size] = static_cast<std::int32_t>(l_filled);
  if (l_filled != l_count) {
    return route_fallback(workspace, sb, cb);
  }

  std::size_t j_filled = 0;
  for (std::size_t ka = 0; ka < lkabtc_size; ++ka) {
    j_offset_host[ka] = static_cast<std::int32_t>(j_filled);
    for (std::size_t j = 0; j < nj_size; ++j) {
      const std::size_t map = kabot0 + ka + j * nkastr_size;
      if (i1[map] != 0) {
        j_index_host[j_filled] = static_cast<std::int32_t>(j);
        j_ja_host[j_filled] = i1[map];
        j_factor_host[j_filled] = xi1s[map];
        ++j_filled;
      }
    }
  }
  j_offset_host[lkabtc_size] = static_cast<std::int32_t>(j_filled);
  if (j_filled != j_count) {
    return route_fallback(workspace, sb, cb);
  }

  double *device_sb = workspace.sb.data;
  const double *device_cb = workspace.cb.data;
  if (shared_match) {
    if (!lucia_sigma_cuda_blocks::sigma_blocks_acquire(sb, cb, sb_count, cb_count,
                                                        &device_sb, &device_cb)) {
      return route_fallback(workspace, sb, cb);
    }
  }
  const bool upload_blocks = local_match && !workspace.session_resident;
  if (!shared_match && (!local_match || upload_blocks) &&
      (!workspace.sb.copy_from(sb, sb_count) || !workspace.cb.copy_from(cb, cb_count))) {
    return route_fallback(workspace, sb, cb);
  }
  if (!shared_match) {
    device_sb = workspace.sb.data;
    device_cb = workspace.cb.data;
  }
  if (!shared_match && (!local_match || upload_blocks) && !workspace.staging.reserve(sb_count)) {
    return route_fallback(workspace, sb, cb);
  }
  if (!slot.upload(payload_bytes)) {
    return route_fallback(workspace, sb, cb);
  }
  workspace.next_route = (workspace.next_route + 1U) % 2U;

  const unsigned char *payload_device = slot.device.data;
  const double *xint_device = reinterpret_cast<const double *>(payload_device + xint_offset);
  const std::int64_t *i3_device = reinterpret_cast<const std::int64_t *>(payload_device + i3_offset);
  const double *xi3s_device = reinterpret_cast<const double *>(payload_device + xi3s_offset);
  const double *inverse_factor_device =
      reinterpret_cast<const double *>(payload_device + inverse_factor_offset);
  const std::int64_t *l_jb_device =
      reinterpret_cast<const std::int64_t *>(payload_device + l_jb_offset);
  const double *l_factor_device =
      reinterpret_cast<const double *>(payload_device + l_factor_offset);
  const std::int64_t *j_ja_device =
      reinterpret_cast<const std::int64_t *>(payload_device + j_ja_offset);
  const double *j_factor_device =
      reinterpret_cast<const double *>(payload_device + j_factor_offset);
  const std::int32_t *inverse_offset_device =
      reinterpret_cast<const std::int32_t *>(payload_device + inverse_offset_offset);
  const std::int32_t *inverse_k_device =
      reinterpret_cast<const std::int32_t *>(payload_device + inverse_k_offset);
  const std::int32_t *inverse_kb_device =
      reinterpret_cast<const std::int32_t *>(payload_device + inverse_kb_offset);
  const std::int32_t *l_offset_device =
      reinterpret_cast<const std::int32_t *>(payload_device + l_offset_offset);
  const std::int32_t *l_index_device =
      reinterpret_cast<const std::int32_t *>(payload_device + l_index_offset);
  const std::int32_t *j_offset_device =
      reinterpret_cast<const std::int32_t *>(payload_device + j_offset_offset);
  const std::int32_t *j_index_device =
      reinterpret_cast<const std::int32_t *>(payload_device + j_index_offset);

  rsbb2bn_cuda_kernel<<<static_cast<unsigned int>(block_count), 256>>>(
      device_sb, device_cb, xint_device,
      i3_device, xi3s_device, inverse_offset_device,
      inverse_k_device, inverse_kb_device, inverse_factor_device, l_offset_device,
      l_index_device, l_jb_device, l_factor_device, j_offset_device,
      j_index_device, j_ja_device, j_factor_device, nib_size, njb_size,
      nkastr_size, kabot0, lkabtc_size, ni_size, nj_size, nk_size, nl_size, ikord == 1);
  bool success = cudaGetLastError() == cudaSuccess;
  if (!shared_match && !local_match && cudaDeviceSynchronize() != cudaSuccess) success = false;

  if (!success) {
    workspace.invalidate_residency();
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
  if (!lucia_cuda::checked_mul(sb_count, sizeof(double), &sb_bytes) ||
      cudaMemcpy(workspace.staging.data, device_sb, sb_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) {
    workspace.invalidate_residency();
    return -1;
  }
  std::memcpy(sb, workspace.staging.data, sb_bytes);
  return 1;
}

extern "C" void lucia_rsbb2bn_cuda_release()
{
  workspace.release();
}
