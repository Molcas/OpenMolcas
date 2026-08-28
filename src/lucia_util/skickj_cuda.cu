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

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>

namespace {

bool checked_mul(std::size_t left, std::size_t right, std::size_t *result) noexcept
{
  if (right != 0 && left > (std::numeric_limits<std::size_t>::max)() / right) {
    return false;
  }
  *result = left * right;
  return true;
}
bool checked_product2(std::size_t first, std::size_t second, std::size_t *result) noexcept
{
  return checked_mul(first, second, result);
}

bool checked_product3(std::size_t first, std::size_t second, std::size_t third, std::size_t *result) noexcept
{
  std::size_t partial = 0;
  return checked_mul(first, second, &partial) && checked_mul(partial, third, result);
}

bool checked_product4(std::size_t first, std::size_t second, std::size_t third, std::size_t fourth, std::size_t *result) noexcept
{
  std::size_t partial = 0;
  return checked_mul(first, second, &partial) && checked_mul(partial, third, &partial) &&
         checked_mul(partial, fourth, result);
}

bool positive_size(int64_t value, std::size_t *result) noexcept
{
  if (value <= 0 || static_cast<std::uintmax_t>(value) >
                    static_cast<std::uintmax_t>((std::numeric_limits<std::size_t>::max)())) {
    return false;
  }
  *result = static_cast<std::size_t>(value);
  return true;
}

template <typename T>
struct DeviceAllocation {
  T *data = nullptr;
  std::size_t capacity = 0;

  bool reserve(std::size_t count) noexcept
  {
    if (count <= capacity) {
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
      data = nullptr;
      capacity = 0;
      cudaFree(replacement);
      return false;
    }
    data = replacement;
    capacity = count;
    return true;
  }

  bool copy_from(const T *host, std::size_t count) noexcept
  {
    if (!reserve(count)) {
      return false;
    }
    std::size_t bytes = 0;
    if (!checked_mul(count, sizeof(T), &bytes)) {
      return false;
    }
    if (cudaMemcpy(data, host, bytes, cudaMemcpyHostToDevice) != cudaSuccess) {
      return false;
    }
    return true;
  }

  bool release() noexcept
  {
    const bool success = data == nullptr || cudaFree(data) == cudaSuccess;
    data = nullptr;
    capacity = 0;
    return success;
  }
};

struct Workspace {
  DeviceAllocation<double> skii;
  DeviceAllocation<double> ckjj;
  DeviceAllocation<double> xijkl;
  DeviceAllocation<int64_t> kbib;
  DeviceAllocation<double> xkbib;
  DeviceAllocation<int64_t> kbjb;
  DeviceAllocation<double> xkbjb;
  double *staging = nullptr;
  std::size_t staging_capacity = 0;
  bool device_ready = false;
  std::size_t max_grid_x = 0;

  bool reserve_staging(std::size_t count) noexcept
  {
    if (count <= staging_capacity) {
      return true;
    }
    std::size_t bytes = 0;
    if (!checked_mul(count, sizeof(double), &bytes)) {
      return false;
    }
    double *replacement = nullptr;
    if (cudaMallocHost(reinterpret_cast<void **>(&replacement), bytes) != cudaSuccess) {
      return false;
    }
    if (staging != nullptr && cudaFreeHost(staging) != cudaSuccess) {
      staging = nullptr;
      staging_capacity = 0;
      cudaFreeHost(replacement);
      return false;
    }
    staging = replacement;
    staging_capacity = count;
    return true;
  }

  bool supports_grid(std::size_t count) noexcept
  {
    if (!device_ready) {
      int device = 0;
      cudaDeviceProp properties{};
      if (cudaGetDevice(&device) != cudaSuccess ||
          cudaGetDeviceProperties(&properties, device) != cudaSuccess || properties.maxGridSize[0] <= 0) {
        return false;
      }
      max_grid_x = static_cast<std::size_t>(properties.maxGridSize[0]);
      device_ready = true;
    }
    return count <= max_grid_x &&
           count <= static_cast<std::size_t>((std::numeric_limits<unsigned int>::max)());
  }

  bool release() noexcept
  {
    bool success = true;
    if (!skii.release()) success = false;
    if (!ckjj.release()) success = false;
    if (!xijkl.release()) success = false;
    if (!kbib.release()) success = false;
    if (!xkbib.release()) success = false;
    if (!kbjb.release()) success = false;
    if (!xkbjb.release()) success = false;
    if (staging != nullptr && cudaFreeHost(staging) != cudaSuccess) success = false;
    staging = nullptr;
    staging_capacity = 0;
    device_ready = false;
    max_grid_x = 0;
    return success;
  }
};

Workspace workspace;

__global__ void skickj_route3_kernel(double *skii, const double *ckjj, const double *xijkl,
                                     int64_t nka, int64_t nkb, int64_t ni, int64_t nj, int64_t nk,
                                     int64_t nl, int64_t maxk, const int64_t *kbib, const double *xkbib,
                                     const int64_t *kbjb, const double *xkbjb, int64_t ikord,
                                     int64_t nib, int64_t njb)
{
  const std::size_t plane = static_cast<std::size_t>(blockIdx.x);
  const std::size_t maxk_size = static_cast<std::size_t>(maxk);
  const std::size_t nka_size = static_cast<std::size_t>(nka);
  const std::size_t ni_size = static_cast<std::size_t>(ni);
  const std::size_t nj_size = static_cast<std::size_t>(nj);

  const int64_t ib = static_cast<int64_t>(plane / ni_size);
  const int64_t i = static_cast<int64_t>(plane % ni_size);
  for (std::size_t ka = static_cast<std::size_t>(threadIdx.x); ka < nka_size; ka += blockDim.x) {
    const std::size_t output_offset = ka + static_cast<std::size_t>(i) * nka_size +
                                      static_cast<std::size_t>(ib) * ni_size * nka_size;
    double result = skii[output_offset];
    for (int64_t kb = 0; kb < nkb; ++kb) {
      for (int64_t k = 0; k < nk; ++k) {
        if (ikord == 1 && i > k) {
          continue;
        }
        const std::size_t map_k_offset = static_cast<std::size_t>(kb) +
                                         static_cast<std::size_t>(k) * maxk_size;
        if (kbib[map_k_offset] != ib + 1) {
          continue;
        }
        const double xkbib_value = xkbib[map_k_offset];
        for (int64_t l = 0; l < nl; ++l) {
          const std::size_t map_l_offset = static_cast<std::size_t>(kb) +
                                           static_cast<std::size_t>(l) * maxk_size;
          const int64_t jb = kbjb[map_l_offset];
          if (jb < 1 || jb > njb) {
            continue;
          }
          const double factor = xkbib_value * xkbjb[map_l_offset];
          double sum = 0.0;
          for (int64_t j = 0; j < nj; ++j) {
            const std::size_t a_offset = ka + static_cast<std::size_t>(j) * nka_size +
                                         static_cast<std::size_t>(jb - 1) * nj_size * nka_size;
            const std::size_t b_offset = static_cast<std::size_t>(j) +
                                         static_cast<std::size_t>(i) * nj_size +
                                         (static_cast<std::size_t>(l) * static_cast<std::size_t>(nk) +
                                          static_cast<std::size_t>(k)) * ni_size * nj_size;
            double b = xijkl[b_offset];
            if (ikord == 1 && i == k) {
              if (j > l && j < nl) {
                b = 0.0;
              } else if (j == l) {
                b *= 0.5;
              }
            }
            sum += ckjj[a_offset] * b;
          }
          result += factor * sum;
        }
      }
    }
    skii[output_offset] = result;
  }
}

} // namespace

extern "C" int64_t lucia_skickj_cuda_route3(
    double *skii, const double *ckjj, const double *xijkl,
    int64_t nka, int64_t nkb, int64_t ni, int64_t nj, int64_t nk, int64_t nl,
    int64_t maxk, const int64_t *kbib, const double *xkbib,
    const int64_t *kbjb, const double *xkbjb, int64_t ikord,
    int64_t nib, int64_t njb, double facs)
{
  if (facs != 1.0 || (ikord != 0 && ikord != 1)) {
    return 0;
  }
  if (skii == nullptr || ckjj == nullptr || xijkl == nullptr || kbib == nullptr ||
      xkbib == nullptr || kbjb == nullptr || xkbjb == nullptr) {
    return 0;
  }
  if (nkb > maxk || (ikord == 1 && (nk > ni || nl > nj))) {
    return 0;
  }

  std::size_t nka_size = 0, nkb_size = 0, ni_size = 0, nj_size = 0, nk_size = 0;
  std::size_t nl_size = 0, maxk_size = 0, nib_size = 0, njb_size = 0;
  if (!positive_size(nka, &nka_size) || !positive_size(nkb, &nkb_size) ||
      !positive_size(ni, &ni_size) || !positive_size(nj, &nj_size) ||
      !positive_size(nk, &nk_size) || !positive_size(nl, &nl_size) ||
      !positive_size(maxk, &maxk_size) || !positive_size(nib, &nib_size) ||
      !positive_size(njb, &njb_size)) {
    return 0;
  }

  std::size_t skii_count = 0, ckjj_count = 0, xijkl_count = 0;
  std::size_t kbib_count = 0, kbjb_count = 0, operation_count = 0, output_plane_count = 0;
  if (!checked_product3(nka_size, ni_size, nib_size, &skii_count) ||
      !checked_product3(nka_size, nj_size, njb_size, &ckjj_count) ||
      !checked_product4(ni_size, nj_size, nk_size, nl_size, &xijkl_count) ||
      !checked_product2(maxk_size, nk_size, &kbib_count) ||
      !checked_product2(maxk_size, nl_size, &kbjb_count) ||
      !checked_product3(nkb_size, nk_size, nl_size, &operation_count) ||
      !checked_product2(nib_size, ni_size, &output_plane_count)) {
    return 0;
  }
  std::size_t work_count = 0;
  constexpr std::size_t minimum_work = 300000000;
  if (!checked_product4(nka_size, ni_size, nj_size, operation_count, &work_count) ||
      work_count < minimum_work) {
    return 0;
  }
  if (operation_count > static_cast<std::size_t>((std::numeric_limits<int64_t>::max)())) {
    return 0;
  }

  if (!workspace.supports_grid(output_plane_count)) {
    return 0;
  }

  std::size_t skii_bytes = 0;
  if (!checked_mul(skii_count, sizeof(double), &skii_bytes)) {
    return 0;
  }

  bool success = true;
  do {
    if (!workspace.skii.copy_from(skii, skii_count) ||
        !workspace.ckjj.copy_from(ckjj, ckjj_count) ||
        !workspace.xijkl.copy_from(xijkl, xijkl_count) ||
        !workspace.kbib.copy_from(kbib, kbib_count) ||
        !workspace.xkbib.copy_from(xkbib, kbib_count) ||
        !workspace.kbjb.copy_from(kbjb, kbjb_count) ||
        !workspace.xkbjb.copy_from(xkbjb, kbjb_count) ||
        !workspace.reserve_staging(skii_count)) {
      success = false;
      break;
    }

    const dim3 block(128, 1, 1);
    const dim3 grid(static_cast<unsigned int>(output_plane_count), 1, 1);
    skickj_route3_kernel<<<grid, block>>>(
        workspace.skii.data, workspace.ckjj.data, workspace.xijkl.data, nka, nkb, ni, nj, nk, nl, maxk,
        workspace.kbib.data, workspace.xkbib.data, workspace.kbjb.data, workspace.xkbjb.data, ikord, nib, njb);
    if (cudaGetLastError() != cudaSuccess) {
      success = false;
    }
    if (cudaDeviceSynchronize() != cudaSuccess) {
      success = false;
    }
    if (!success) {
      break;
    }

    if (cudaMemcpy(workspace.staging, workspace.skii.data, skii_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) {
      success = false;
      break;
    }
  } while (false);

  if (success) {
    std::memcpy(skii, workspace.staging, skii_bytes);
  }
  return success ? 1 : 0;
}

extern "C" void lucia_skickj_cuda_release()
{
  workspace.release();
}
