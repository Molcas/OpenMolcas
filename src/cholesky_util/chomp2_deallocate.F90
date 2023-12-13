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

use ChoMP2, only: ChoMP2_allocated, iFirst, iFirstS, LiMatij, LiPQprod, LiT1am, LnBatOrb, LnMatij, LnOcc, LnPQprod, LnT1am, lUnit, &
                  NumBatOrb, NumOcc
use stdalloc, only: mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: irc

irc = 0

call ChoMP2g_deallocate(irc)

if (ChoMP2_allocated) then

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

end if

end subroutine ChoMP2_deallocate
