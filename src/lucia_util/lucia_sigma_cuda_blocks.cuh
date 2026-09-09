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


#ifndef LUCIA_SIGMA_CUDA_BLOCKS_CUH
#define LUCIA_SIGMA_CUDA_BLOCKS_CUH

#include <cstddef>

namespace lucia_sigma_cuda_blocks {

bool sigma_blocks_match(double *sb, const double *cb, std::size_t sb_count,
                        std::size_t cb_count) noexcept;

bool sigma_blocks_acquire(double *sb, const double *cb, std::size_t sb_count,
                          std::size_t cb_count, double **sb_device,
                          const double **cb_device) noexcept;

bool sigma_blocks_flush_if_match(double *sb, const double *cb,
                                 bool *matched) noexcept;

bool sigma_blocks_flush_active(bool *active) noexcept;

} // namespace lucia_sigma_cuda_blocks

#endif
