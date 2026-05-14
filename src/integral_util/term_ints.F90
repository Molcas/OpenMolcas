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
! Copyright (C) 1998, Roland Lindh                                     *
!***********************************************************************

subroutine Term_Ints()
!***********************************************************************
!                                                                      *
!     Object: to deallocate memory in association with two-electron    *
!             calculations.                                            *
!                                                                      *
!     Author: Roland Lindh, Chemical Physics, University of Lund,      *
!             Sweden. January '98.                                     *
!***********************************************************************

use k2_arrays, only: Aux, Destroy_BraKet_Base, FT, iSOSym
use iSD_Data, only: iCntr, iCntr, iSh2Sh, iShOff, iSO2Sh, nShBf
use stdalloc, only: mma_deallocate

implicit none

!                                                                      *
!***********************************************************************
!                                                                      *
! In case of semi-direct mode the memory is released externally.

call RlsMem_Ints()
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(FT,safe='*')

call Destroy_Braket_base()

call mma_deallocate(Aux,safe='*')

call mma_deallocate(iSOSym,safe='*')
!                                                                      *
!***********************************************************************
!                                                                      *
if (allocated(nShBf)) then
  call mma_deallocate(nShBF)
  call mma_deallocate(iShOff)
  call mma_deallocate(iSh2Sh)
  call mma_deallocate(iSO2Sh)
  call mma_deallocate(iCntr)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Free memory for K2 data

call FreeK2()
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Term_Ints
