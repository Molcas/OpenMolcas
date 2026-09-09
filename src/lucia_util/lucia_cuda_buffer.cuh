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


#ifndef LUCIA_CUDA_BUFFER_CUH
#define LUCIA_CUDA_BUFFER_CUH

#include <cuda_runtime.h>

#include <cstddef>
#include <cstdint>
#include <limits>

namespace lucia_cuda {

inline bool checked_mul(std::size_t left, std::size_t right, std::size_t *result) noexcept
{
  if (right != 0 && left > (std::numeric_limits<std::size_t>::max)() / right) {
    return false;
  }
  *result = left * right;
  return true;
}

inline bool positive_size(std::int64_t value, std::size_t *result) noexcept
{
  if (value <= 0 || static_cast<std::uintmax_t>(value) >
                    static_cast<std::uintmax_t>((std::numeric_limits<std::size_t>::max)())) {
    return false;
  }
  *result = static_cast<std::size_t>(value);
  return true;
}

template <typename T>
struct DeviceBuffer {
  T *data = nullptr;
  std::size_t capacity = 0;

  DeviceBuffer() = default;
  DeviceBuffer(const DeviceBuffer &) = delete;
  DeviceBuffer &operator=(const DeviceBuffer &) = delete;

  bool reserve(std::size_t count) noexcept
  {
    if (count <= capacity && data != nullptr) {
      return true;
    }
    if (count == 0) {
      return true;
    }

    std::size_t bytes = 0;
    if (!checked_mul(count, sizeof(T), &bytes)) {
      return false;
    }
    T *replacement = nullptr;
    if (cudaMalloc(reinterpret_cast<void **>(&replacement), bytes) != cudaSuccess) {
      return false;
    }
    if (data != nullptr && cudaFree(data) != cudaSuccess) {
      cudaFree(replacement);
      return false;
    }
    data = replacement;
    capacity = count;
    return true;
  }

  bool copy_from(const T *host, std::size_t count) noexcept
  {
    if (count == 0) {
      return true;
    }
    if (host == nullptr || !reserve(count)) {
      return false;
    }
    std::size_t bytes = 0;
    if (!checked_mul(count, sizeof(T), &bytes)) {
      return false;
    }
    return cudaMemcpy(data, host, bytes, cudaMemcpyHostToDevice) == cudaSuccess;
  }

  bool release() noexcept
  {
    const bool success = data == nullptr || cudaFree(data) == cudaSuccess;
    data = nullptr;
    capacity = 0;
    return success;
  }
};

template <typename T>
struct PinnedBuffer {
  T *data = nullptr;
  std::size_t capacity = 0;

  PinnedBuffer() = default;
  PinnedBuffer(const PinnedBuffer &) = delete;
  PinnedBuffer &operator=(const PinnedBuffer &) = delete;

  bool reserve(std::size_t count) noexcept
  {
    if (count <= capacity && data != nullptr) {
      return true;
    }
    if (count == 0) {
      return true;
    }

    std::size_t bytes = 0;
    if (!checked_mul(count, sizeof(T), &bytes)) {
      return false;
    }
    T *replacement = nullptr;
    if (cudaMallocHost(reinterpret_cast<void **>(&replacement), bytes) != cudaSuccess) {
      return false;
    }
    if (data != nullptr && cudaFreeHost(data) != cudaSuccess) {
      cudaFreeHost(replacement);
      return false;
    }
    data = replacement;
    capacity = count;
    return true;
  }

  bool release() noexcept
  {
    const bool success = data == nullptr || cudaFreeHost(data) == cudaSuccess;
    data = nullptr;
    capacity = 0;
    return success;
  }
};

} // namespace lucia_cuda

#endif
