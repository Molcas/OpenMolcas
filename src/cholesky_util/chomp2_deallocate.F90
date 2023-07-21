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
! Copyright (C) 2004,2005, Thomas Bondo Pedersen                       *
!               2010, Jonas Bostrom                                    *
!               2021, Roland Lindh                                     *
!***********************************************************************

subroutine ChoMP2_deallocate(irc)
!
! Purpose: to deallocate memory of the  Cholesky MP2 program.

use ChoMP2, only: ChoMP2_allocated
use ChoMP2, only: iFirst, iFirstS, NumOcc, LnOcc, LnT1am, LiT1am
use ChoMP2, only: LnMatij, LiMatij, lUnit, NumBatOrb, LnBatOrb
use ChoMP2, only: LnPQprod, LiPQprod
use stdalloc

implicit real*8(a-h,o-z)
#include "chomp2.fh"

irc = 0

call ChoMP2g_deallocate(irc)

if (.not. ChoMP2_allocated) return

call mma_deallocate(LiPQprod)
call mma_deallocate(LnPQprod)
call mma_deallocate(LnBatOrb)
call mma_deallocate(NumBatOrb)
call mma_deallocate(lUnit)
call mma_deallocate(LiMatij)
call mma_deallocate(LnMatij)
call mma_deallocate(LiT1am)
call mma_deallocate(LnT1am)
call mma_deallocate(LnOcc)
call mma_deallocate(NumOcc)
call mma_deallocate(iFirstS)
call mma_deallocate(iFirst)

ChoMP2_allocated = .false.

end subroutine ChoMP2_deallocate
