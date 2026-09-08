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

struct HostRegion {
  const void *base = nullptr;
  std::size_t bytes = 0;
  bool registered = false;
  bool owned = false;

  void clear() noexcept {
    base = nullptr;
    bytes = 0;
    registered = false;
    owned = false;
  }
};

bool range_contains(const HostRegion &region, const void *pointer, std::size_t bytes) noexcept {
  if (!region.registered || region.base == nullptr || pointer == nullptr) {
    return false;
  }
  const std::uintptr_t base = reinterpret_cast<std::uintptr_t>(region.base);
  const std::uintptr_t address = reinterpret_cast<std::uintptr_t>(pointer);
  const std::uintptr_t maximum = (std::numeric_limits<std::uintptr_t>::max)();
  if (base > maximum - region.bytes || address > maximum - bytes) {
    return false;
  }
  const std::uintptr_t region_end = base + region.bytes;
  const std::uintptr_t address_end = address + bytes;
  return address >= base && address_end <= region_end;
}

bool register_host_region(HostRegion *region, void *base, std::size_t bytes) noexcept {
  if (region == nullptr || base == nullptr || bytes == 0) {
    return false;
  }
  if (cudaHostRegister(base, bytes, cudaHostRegisterPortable) != cudaSuccess) {
    return false;
  }
  region->base = base;
  region->bytes = bytes;
  region->registered = true;
  region->owned = true;
  return true;
}

bool unregister_host_region(HostRegion *region) noexcept {
  if (region == nullptr || !region->registered) {
    if (region != nullptr) {
      region->clear();
    }
    return true;
  }
  if (!region->owned) {
    region->clear();
    return true;
  }
  if (cudaHostUnregister(const_cast<void *>(region->base)) != cudaSuccess) {
    return false;
  }
  region->clear();
  return true;
}

struct Manager {
  lucia_cuda::DeviceBuffer<double> sb;
  lucia_cuda::DeviceBuffer<double> cb;
  lucia_cuda::PinnedBuffer<double> staging;
  bool active = false;
  bool resident = false;
  bool zero_sb = false;
  double *host_sb = nullptr;
  const double *host_cb = nullptr;
  std::size_t sb_count = 0;
  std::size_t cb_count = 0;
  HostRegion host_sb_region;
  HostRegion host_cb_region;
  bool host_active = false;

  bool matches(double *sb_host, const double *cb_host, std::size_t requested_sb, std::size_t requested_cb) const noexcept {
    return active && host_sb == sb_host && host_cb == cb_host && requested_sb > 0 && requested_cb > 0 && requested_sb <= sb_count
           && requested_cb <= cb_count;
  }

  bool pointer_match(double *sb_host, const double *cb_host) const noexcept {
    return active && host_sb == sb_host && host_cb == cb_host;
  }

  bool begin(double *sb_host, const double *cb_host, std::size_t full_sb, std::size_t full_cb,
             bool zero_sb_device = false) noexcept {
    if (active) {
      return false;
    }
    active = true;
    resident = false;
    zero_sb = zero_sb_device;
    host_sb = sb_host;
    host_cb = cb_host;
    sb_count = full_sb;
    cb_count = full_cb;
    return true;
  }

  bool acquire(double *sb_host, const double *cb_host, std::size_t requested_sb, std::size_t requested_cb, double **sb_device,
               const double **cb_device) noexcept {
    if (sb_device == nullptr || cb_device == nullptr || !matches(sb_host, cb_host, requested_sb, requested_cb)) {
      return false;
    }
    if (!resident) {
      std::size_t sb_bytes = 0;
      const bool sb_ready = zero_sb ? lucia_cuda::checked_mul(sb_count, sizeof(double), &sb_bytes) && sb.reserve(sb_count)
                                          && cudaMemset(sb.data, 0, sb_bytes) == cudaSuccess
                                    : sb.copy_from(host_sb, sb_count);
      if (!sb_ready || !cb.copy_from(host_cb, cb_count)) {
        resident = false;
        return false;
      }
      zero_sb = false;
      resident = true;
    }
    *sb_device = sb.data;
    *cb_device = cb.data;
    return true;
  }

  bool flush() noexcept {
    if (!active) {
      return true;
    }
    if (!resident) {
      zero_sb = false;
      return true;
    }
    resident = false;
    std::size_t bytes = 0;
    if (!lucia_cuda::checked_mul(sb_count, sizeof(double), &bytes)) {
      return false;
    }
    if (range_contains(host_sb_region, host_sb, bytes)) {
      return cudaMemcpy(host_sb, sb.data, bytes, cudaMemcpyDeviceToHost) == cudaSuccess;
    }
    if (!staging.reserve(sb_count) || cudaMemcpy(staging.data, sb.data, bytes, cudaMemcpyDeviceToHost) != cudaSuccess) {
      return false;
    }
    std::memcpy(host_sb, staging.data, bytes);
    return true;
  }

  std::int64_t host_begin(double *sb_host, const double *cb_host, std::int64_t requested_sb, std::int64_t requested_cb) noexcept {
    if (host_active || sb_host == nullptr || cb_host == nullptr) {
      return 0;
    }
    std::size_t sb_size = 0;
    std::size_t cb_size = 0;
    std::size_t sb_bytes = 0;
    std::size_t cb_bytes = 0;
    if (!lucia_cuda::positive_size(requested_sb, &sb_size) || !lucia_cuda::positive_size(requested_cb, &cb_size)
        || !lucia_cuda::checked_mul(sb_size, sizeof(double), &sb_bytes)
        || !lucia_cuda::checked_mul(cb_size, sizeof(double), &cb_bytes)) {
      return 0;
    }

    HostRegion sb_region;
    HostRegion cb_region;
    if (!register_host_region(&sb_region, sb_host, sb_bytes)) {
      return 0;
    }
    if (!register_host_region(&cb_region, const_cast<double *>(cb_host), cb_bytes)) {
      if (!unregister_host_region(&sb_region)) {
        host_sb_region = sb_region;
        host_cb_region.clear();
        host_active = true;
        return -1;
      }
      return 0;
    }

    host_sb_region = sb_region;
    host_cb_region = cb_region;
    host_active = true;
    return 1;
  }

  bool host_end() noexcept {
    if (!host_active) {
      return true;
    }
    if (active) {
      return false;
    }
    bool success = unregister_host_region(&host_sb_region);
    if (!unregister_host_region(&host_cb_region)) {
      success = false;
    }
    if (success) {
      host_active = false;
      host_sb_region.clear();
      host_cb_region.clear();
    }
    return success;
  }

  void clear() noexcept {
    active = false;
    resident = false;
    zero_sb = false;
    host_sb = nullptr;
    host_cb = nullptr;
    sb_count = 0;
    cb_count = 0;
  }

  bool end() noexcept {
    const bool success = flush();
    clear();
    return success;
  }

  bool release() noexcept {
    bool success = end();
    if (!host_end())
      success = false;
    if (!sb.release())
      success = false;
    if (!cb.release())
      success = false;
    if (!staging.release())
      success = false;
    return success;
  }
};

Manager manager;

} // namespace

namespace lucia_sigma_cuda_blocks {

bool sigma_blocks_match(double *sb, const double *cb, std::size_t sb_count, std::size_t cb_count) noexcept {
  return manager.matches(sb, cb, sb_count, cb_count);
}

bool sigma_blocks_acquire(double *sb, const double *cb, std::size_t sb_count, std::size_t cb_count, double **sb_device,
                          const double **cb_device) noexcept {
  return manager.acquire(sb, cb, sb_count, cb_count, sb_device, cb_device);
}

bool sigma_blocks_flush_if_match(double *sb, const double *cb, bool *matched) noexcept {
  if (matched == nullptr) {
    return false;
  }
  *matched = manager.pointer_match(sb, cb);
  return !*matched || manager.flush();
}

bool sigma_blocks_flush_active(bool *active) noexcept {
  if (active == nullptr) {
    return false;
  }
  *active = manager.active;
  return !*active || manager.flush();
}

} // namespace lucia_sigma_cuda_blocks

extern "C" std::int64_t lucia_sigma_cuda_blocks_begin(double *sb, const double *cb, std::int64_t sb_count, std::int64_t cb_count) {
  if (sb == nullptr || cb == nullptr || manager.active) {
    return 0;
  }
  std::size_t sb_size = 0;
  std::size_t cb_size = 0;
  if (!lucia_cuda::positive_size(sb_count, &sb_size) || !lucia_cuda::positive_size(cb_count, &cb_size)
      || !manager.begin(sb, cb, sb_size, cb_size)) {
    return 0;
  }
  return 1;
}

extern "C" std::int64_t lucia_sigma_cuda_blocks_begin_zeroed(double *sb, const double *cb, std::int64_t sb_count,
                                                             std::int64_t cb_count) {
  if (sb == nullptr || cb == nullptr || manager.active) {
    return 0;
  }
  std::size_t sb_size = 0;
  std::size_t cb_size = 0;
  if (!lucia_cuda::positive_size(sb_count, &sb_size) || !lucia_cuda::positive_size(cb_count, &cb_size)
      || !manager.begin(sb, cb, sb_size, cb_size, true)) {
    return 0;
  }
  return 1;
}

extern "C" std::int64_t lucia_sigma_cuda_blocks_end() {
  return manager.end() ? 1 : -1;
}

extern "C" std::int64_t lucia_sigma_cuda_blocks_host_begin(double *sb, const double *cb, std::int64_t sb_count,
                                                           std::int64_t cb_count) {
  return manager.host_begin(sb, cb, sb_count, cb_count);
}

extern "C" std::int64_t lucia_sigma_cuda_blocks_host_end() {
  return manager.host_end() ? 1 : -1;
}

extern "C" void lucia_sigma_cuda_blocks_release() {
  manager.release();
}
