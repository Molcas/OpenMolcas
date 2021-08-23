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

subroutine ProdsS_1(AB,iAB,CMO,nMO,Y)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iAB, nMO
real(kind=wp), intent(in) :: AB(iAB*(iAB+1)/2), CMO(iAB,nMO)
real(kind=wp), intent(out) :: Y(iAB,nMO)
real(kind=wp), allocatable :: ABSq(:,:)

call mma_allocate(ABSq,iAB,iAB,Label='ABSq')
call SQUARE(AB,ABSq,1,iAB,iAB)
call DGEMM_('N','N',iAB,nMO,iAB,One,ABSq,iAB,CMO,iAB,Zero,Y,iAB)
call mma_deallocate(ABSq)

return

end subroutine ProdsS_1
