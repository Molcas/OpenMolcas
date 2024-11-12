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

subroutine Free_HerRW()

use Her_RW, only: HerR, HerW, iHerR, iHerw
use stdalloc, only: mma_deallocate

implicit none

call mma_deallocate(iHerR,safe='*')
call mma_deallocate(iHerW,safe='*')
call mma_deallocate(HerR,safe='*')
call mma_deallocate(HerW,safe='*')

return

end subroutine Free_HerRW
