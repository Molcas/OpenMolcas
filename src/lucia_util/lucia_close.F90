
!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Lucia_Close()
use definitions, only: iwp

use lucia_data, only: LUC, LUDIA, LUHC, LUMOUT, LUSC1, LUSC2, LUSC3, LUSC34, LUSC35, LUSC36, LUSC37, LUSC38, LUSC39, LUSC40

implicit none
logical(kind=iwp), external :: is_opened

! Free memory allocated by Lucia

call FREESTR_GAS()
#ifdef _CUDA_BLAS_
CudaStatus = LUCIA_SIGMA_CUDA_BLOCKS_HOST_END()
if (CudaStatus == -1_c_int64_t) call SYSABENDMSG('lucia_util/lucia_close','CUDA execution failed','')
#endif
call DeAlloc_Lucia()

! Close any files opened by Lucia

If (is_opened(LUDIA )) call DAClos(LUDIA)
If (is_opened(LUC   )) call DAClos(LUC)
If (is_opened(LUHC  )) call DAClos(LUHC)
If (is_opened(LUSC1 )) call DAClos(LUSC1)
If (is_opened(LUSC2 )) call DAClos(LUSC2)
If (is_opened(LUSC3 )) call DAClos(LUSC3)
If (is_opened(LUSC34)) call DAClos(LUSC34)
If (is_opened(LUSC35)) call DAClos(LUSC35)
If (is_opened(LUSC36)) call DAClos(LUSC36)
If (is_opened(LUSC37)) call DAClos(LUSC37)
If (is_opened(LUSC38)) call DAClos(LUSC38)
If (is_opened(LUSC39)) call DAClos(LUSC39)
If (is_opened(LUSC40)) call DAClos(LUSC40)
If (is_opened(LUMOUT)) call DAClos(LUMOUT)

end subroutine Lucia_Close
