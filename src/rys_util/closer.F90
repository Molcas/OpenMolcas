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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************

subroutine CloseR()
!***********************************************************************
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             September '90                                            *
!***********************************************************************

use vRys_RW, only: Cff, ddx, HerR2, HerW2, iHerR2, iHerW2, Map, TMax, x0
use abdata, only: close_abdata
#ifdef _RYS_SCRATCH_
use RysScratch, only: UnSetAux
#endif
use stdalloc, only: mma_deallocate

implicit none

#ifdef _RYS_SCRATCH_
call UnSetAux()
#endif

call close_abdata()

if (.not. allocated(iHerW2)) return
call mma_deallocate(iHerW2)
call mma_deallocate(iHerR2)
call mma_deallocate(HerW2)
call mma_deallocate(HerR2)
call mma_deallocate(Cff)
call mma_deallocate(x0)
call mma_deallocate(Map)
call mma_deallocate(ddx)
call mma_deallocate(TMax)

return

end subroutine CloseR
