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

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>

namespace {

constexpr std::size_t tile_size = 16;
constexpr std::size_t maps_per_split = 16;
constexpr std::size_t max_split_count = 16;
constexpr std::size_t minimum_active_work = 10000000;

bool checked_add(std::size_t left, std::size_t right, std::size_t *result) noexcept
{
  if (left > (std::numeric_limits<std::size_t>::max)() - right) {
    return false;
  }
  *result = left + right;
  return true;
}
struct Workspace {
  lucia_cuda::DeviceBuffer<double> rho2b;
  lucia_cuda::DeviceBuffer<double> skii;
  lucia_cuda::DeviceBuffer<double> ckjj;
  lucia_cuda::DeviceBuffer<std::int64_t> kbib;
  lucia_cuda::DeviceBuffer<double> xkbib;
  lucia_cuda::DeviceBuffer<std::int64_t> kbjb;
  lucia_cuda::DeviceBuffer<double> xkbjb;
  lucia_cuda::PinnedBuffer<double> staging;
  bool device_ready = false;
  std::size_t max_grid_x = 0;
  std::size_t max_grid_y = 0;
  std::size_t max_grid_z = 0;

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

  bool release() noexcept
  {
    bool success = true;
    if (!rho2b.release()) success = false;
    if (!skii.release()) success = false;
    if (!ckjj.release()) success = false;
    if (!kbib.release()) success = false;
    if (!xkbib.release()) success = false;
    if (!kbjb.release()) success = false;
    if (!xkbjb.release()) success = false;
    if (!staging.release()) success = false;
    device_ready = false;
    max_grid_x = 0;
    max_grid_y = 0;
    max_grid_z = 0;
    return success;
  }
};

Workspace workspace;

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

__global__ void abtor2_cuda_kernel(double *rho2b, const double *skii, const double *ckjj,
                                   const std::int64_t *kbib, const double *xkbib,
                                   const std::int64_t *kbjb, const double *xkbjb,
                                   std::int64_t nka, std::int64_t nkb, std::int64_t ni,
                                   std::int64_t nj, std::int64_t nk, std::int64_t nl,
                                   std::int64_t maxk, std::int64_t nib, std::int64_t njb,
                                   std::size_t split_count, std::size_t maps_in_split)
{
  __shared__ double a_tile[tile_size][tile_size + 1];
  __shared__ double b_tile[tile_size][tile_size + 1];
  __shared__ std::int64_t ib;
  __shared__ std::int64_t jb;
  __shared__ double factor;

  const std::size_t nka_size = static_cast<std::size_t>(nka);
  const std::size_t ni_size = static_cast<std::size_t>(ni);
  const std::size_t nj_size = static_cast<std::size_t>(nj);
  const std::size_t nk_size = static_cast<std::size_t>(nk);
  const std::size_t maxk_size = static_cast<std::size_t>(maxk);
  const std::size_t plane_split = static_cast<std::size_t>(blockIdx.z);
  const std::size_t plane = plane_split / split_count;
  const std::size_t split = plane_split - plane * split_count;
  const std::size_t k = plane % nk_size;
  const std::size_t l = plane / nk_size;
  const std::size_t i = static_cast<std::size_t>(blockIdx.x) * tile_size + threadIdx.y;
  const std::size_t j = static_cast<std::size_t>(blockIdx.y) * tile_size + threadIdx.x;
  const std::size_t map_begin = split * maps_in_split;
  const std::size_t map_end = map_begin + maps_in_split < static_cast<std::size_t>(nkb) ?
                              map_begin + maps_in_split : static_cast<std::size_t>(nkb);
  double result = 0.0;
  for (std::size_t kb = map_begin; kb < map_end; ++kb) {
    const std::size_t k_offset = kb + k * maxk_size;
    const std::size_t l_offset = kb + l * maxk_size;
    if (threadIdx.x == 0 && threadIdx.y == 0) {
      ib = kbib[k_offset];
      jb = kbjb[l_offset];
      if (ib < 1 || ib > nib || jb < 1 || jb > njb) {
        ib = 0;
        jb = 0;
      } else {
        factor = xkbib[k_offset] * xkbjb[l_offset];
      }
    }
    __syncthreads();
    if (ib != 0 && jb != 0) {
      double product = 0.0;
      for (std::size_t base = 0; base < nka_size; base += tile_size) {
        const std::size_t a_ka = base + threadIdx.x;
        const std::size_t b_ka = base + threadIdx.y;
        a_tile[threadIdx.y][threadIdx.x] = i < ni_size && a_ka < nka_size ?
            skii[a_ka + i * nka_size + static_cast<std::size_t>(ib - 1) * ni_size * nka_size] : 0.0;
        b_tile[threadIdx.y][threadIdx.x] = j < nj_size && b_ka < nka_size ?
            ckjj[b_ka + j * nka_size + static_cast<std::size_t>(jb - 1) * nj_size * nka_size] : 0.0;
        __syncthreads();
        if (i < ni_size && j < nj_size) {
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

  if (i < ni_size && j < nj_size && result != 0.0) {
    const std::size_t output = i + j * ni_size + (k + l * nk_size) * ni_size * nj_size;
    atomic_add(rho2b + output, result);
  }
}

} // namespace

extern "C" int64_t lucia_abtor2_cuda_route(
    double *rho2b, const double *skii, const double *ckjj,
    int64_t nka, int64_t nkb, int64_t ni, int64_t nj, int64_t nk, int64_t nl,
    int64_t maxk, const int64_t *kbib, const double *xkbib,
    const int64_t *kbjb, const double *xkbjb, int64_t ikord,
    int64_t nib, int64_t njb)
{
  if (ikord != 0 || rho2b == nullptr || skii == nullptr || ckjj == nullptr ||
      kbib == nullptr || xkbib == nullptr || kbjb == nullptr || xkbjb == nullptr) {
    return 0;
  }

  std::size_t nka_size = 0;
  std::size_t nkb_size = 0;
  std::size_t ni_size = 0;
  std::size_t nj_size = 0;
  std::size_t nk_size = 0;
  std::size_t nl_size = 0;
  std::size_t maxk_size = 0;
  std::size_t nib_size = 0;
  std::size_t njb_size = 0;
  if (!lucia_cuda::positive_size(nka, &nka_size) ||
      !lucia_cuda::positive_size(nkb, &nkb_size) ||
      !lucia_cuda::positive_size(ni, &ni_size) ||
      !lucia_cuda::positive_size(nj, &nj_size) ||
      !lucia_cuda::positive_size(nk, &nk_size) ||
      !lucia_cuda::positive_size(nl, &nl_size) ||
      !lucia_cuda::positive_size(maxk, &maxk_size) ||
      !lucia_cuda::positive_size(nib, &nib_size) ||
      !lucia_cuda::positive_size(njb, &njb_size) || nkb_size > maxk_size) {
    return 0;
  }

  std::size_t kbib_count = 0;
  if (!lucia_cuda::checked_mul(maxk_size, nk_size, &kbib_count)) {
    return 0;
  }
  std::size_t kbjb_count = 0;
  if (!lucia_cuda::checked_mul(maxk_size, nl_size, &kbjb_count)) {
    return 0;
  }

  // Route only sufficiently large active ABTOR2 work to CUDA.
  std::size_t active_pairs = 0;
  for (std::size_t kb = 0; kb < nkb_size; ++kb) {
    std::size_t active_k = 0;
    for (std::size_t k = 0; k < nk_size; ++k) {
      if (kbib[kb + k * maxk_size] != 0) ++active_k;
    }
    std::size_t active_l = 0;
    for (std::size_t l = 0; l < nl_size; ++l) {
      if (kbjb[kb + l * maxk_size] != 0) ++active_l;
    }
    std::size_t kb_pairs = 0;
    if (!lucia_cuda::checked_mul(active_k, active_l, &kb_pairs) ||
        !checked_add(active_pairs, kb_pairs, &active_pairs)) {
      return 0;
    }
  }
  std::size_t active_work = active_pairs;
  if (!lucia_cuda::checked_mul(active_work, ni_size, &active_work) ||
      !lucia_cuda::checked_mul(active_work, nj_size, &active_work) ||
      !lucia_cuda::checked_mul(active_work, nka_size, &active_work) ||
      active_work < minimum_active_work) {
    return 0;
  }

  std::size_t partial = 0;
  std::size_t rho2b_count = 0;
  if (!lucia_cuda::checked_mul(ni_size, nj_size, &partial) ||
      !lucia_cuda::checked_mul(partial, nk_size, &partial) ||
      !lucia_cuda::checked_mul(partial, nl_size, &rho2b_count)) {
    return 0;
  }

  std::size_t skii_count = 0;
  if (!lucia_cuda::checked_mul(nka_size, ni_size, &partial) ||
      !lucia_cuda::checked_mul(partial, nib_size, &skii_count)) {
    return 0;
  }

  std::size_t ckjj_count = 0;
  if (!lucia_cuda::checked_mul(nka_size, nj_size, &partial) ||
      !lucia_cuda::checked_mul(partial, njb_size, &ckjj_count)) {
    return 0;
  }

  std::size_t grid_x = ni_size / tile_size + (ni_size % tile_size != 0);
  std::size_t grid_y = nj_size / tile_size + (nj_size % tile_size != 0);
  std::size_t split_count = nkb_size / maps_per_split + (nkb_size % maps_per_split != 0);
  if (split_count > max_split_count) {
    split_count = max_split_count;
  }
  const std::size_t maps_in_split = nkb_size / split_count + (nkb_size % split_count != 0);
  std::size_t plane_count = 0;
  std::size_t grid_z = 0;
  if (!lucia_cuda::checked_mul(nk_size, nl_size, &plane_count) ||
      !lucia_cuda::checked_mul(plane_count, split_count, &grid_z) ||
      !workspace.supports_grid(grid_x, grid_y, grid_z)) {
    return 0;
  }

  std::size_t rho2b_bytes = 0;
  if (!lucia_cuda::checked_mul(rho2b_count, sizeof(double), &rho2b_bytes)) {
    return 0;
  }

  if (!workspace.rho2b.copy_from(rho2b, rho2b_count) ||
      !workspace.skii.copy_from(skii, skii_count) ||
      !workspace.ckjj.copy_from(ckjj, ckjj_count) ||
      !workspace.kbib.copy_from(kbib, kbib_count) ||
      !workspace.xkbib.copy_from(xkbib, kbib_count) ||
      !workspace.kbjb.copy_from(kbjb, kbjb_count) ||
      !workspace.xkbjb.copy_from(xkbjb, kbjb_count) ||
      !workspace.staging.reserve(rho2b_count)) {
    return 0;
  }

  const dim3 block(static_cast<unsigned int>(tile_size), static_cast<unsigned int>(tile_size), 1);
  const dim3 grid(static_cast<unsigned int>(grid_x), static_cast<unsigned int>(grid_y),
                  static_cast<unsigned int>(grid_z));
  abtor2_cuda_kernel<<<grid, block>>>(
      workspace.rho2b.data, workspace.skii.data, workspace.ckjj.data,
      workspace.kbib.data, workspace.xkbib.data, workspace.kbjb.data, workspace.xkbjb.data,
      nka, nkb, ni, nj, nk, nl, maxk, nib, njb, split_count, maps_in_split);
  bool success = true;
  if (cudaGetLastError() != cudaSuccess) {
    success = false;
  }
  if (cudaDeviceSynchronize() != cudaSuccess) {
    success = false;
  }
  if (!success) {
    return 0;
  }
  if (cudaMemcpy(workspace.staging.data, workspace.rho2b.data, rho2b_bytes,
                 cudaMemcpyDeviceToHost) != cudaSuccess) {
    return 0;
  }

  std::memcpy(rho2b, workspace.staging.data, rho2b_bytes);
  return 1;
}

extern "C" void lucia_abtor2_cuda_release()
{
  workspace.release();
}
