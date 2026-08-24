!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2026, Meng Wang                                        *
!***********************************************************************

subroutine Lucia_Close()

use lucia_data, only: LUC, LUDIA, LUHC, LUMOUT, LUSC1, LUSC2, LUSC3, LUSC34, LUSC35, LUSC36, LUSC37, LUSC38, LUSC39, LUSC40
#ifdef _CUDA_BLAS_
use SKICKJ_CUDA_INTERFACE, only: LUCIA_SKICKJ_CUDA_RELEASE
use ABTOR2_CUDA_INTERFACE, only: LUCIA_ABTOR2_CUDA_RELEASE
use RSBB1E_CUDA_INTERFACE, only: LUCIA_RSBB1E_CUDA_RELEASE
use LUCIA_SIGMA_CUDA_BLOCKS_INTERFACE, only: LUCIA_SIGMA_CUDA_BLOCKS_HOST_END, &
                                             LUCIA_SIGMA_CUDA_BLOCKS_RELEASE
use GSBBD1_CUDA_INTERFACE, only: LUCIA_GSBBD1_CUDA_RELEASE
use GSBBD2A_CUDA_INTERFACE, only: LUCIA_GSBBD2A_CUDA_RELEASE
use RSBB2A_CUDA_INTERFACE, only: LUCIA_RSBB2A_CUDA_RELEASE
use RSBB2BN_CUDA_INTERFACE, only: LUCIA_RSBB2BN_CUDA_RELEASE
use GSBBD2B_CUDA_INTERFACE, only: LUCIA_GSBBD2B_CUDA_RELEASE
use, intrinsic :: iso_c_binding, only: c_int64_t
#endif

implicit none
#ifdef _CUDA_BLAS_
integer(c_int64_t) :: CudaStatus
#endif

! Free memory allocated by Lucia

call FREESTR_GAS()
#ifdef _CUDA_BLAS_
CudaStatus = LUCIA_SIGMA_CUDA_BLOCKS_HOST_END()
if (CudaStatus == -1_c_int64_t) call SYSABENDMSG('lucia_util/lucia_close','CUDA execution failed','')
#endif
call DeAlloc_Lucia()

! Close any files opened by Lucia

call DAClos(LUDIA)
call DAClos(LUC)
call DAClos(LUHC)
call DAClos(LUSC1)
call DAClos(LUSC2)
call DAClos(LUSC3)
call DAClos(LUSC34)
call DAClos(LUSC35)
call DAClos(LUSC36)
call DAClos(LUSC37)
call DAClos(LUSC38)
call DAClos(LUSC39)
call DAClos(LUSC40)
call DAClos(LUMOUT)
#ifdef _CUDA_BLAS_
call LUCIA_SKICKJ_CUDA_RELEASE()
call LUCIA_ABTOR2_CUDA_RELEASE()
call LUCIA_RSBB1E_CUDA_RELEASE()
call LUCIA_GSBBD1_CUDA_RELEASE()
call LUCIA_GSBBD2A_CUDA_RELEASE()
call LUCIA_RSBB2A_CUDA_RELEASE()
call LUCIA_SIGMA_CUDA_BLOCKS_RELEASE()
call LUCIA_RSBB2BN_CUDA_RELEASE()
call LUCIA_GSBBD2B_CUDA_RELEASE()
#endif

end subroutine Lucia_Close
