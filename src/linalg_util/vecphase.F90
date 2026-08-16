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

subroutine VecPhase(A,nA)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nA
real(kind=wp), intent(inout) :: A(nA)
integer(kind=iwp) :: i
real(kind=wp) :: Rand
real(kind=wp), external :: Random_Molcas
integer(kind=iwp) :: iSeed=17

! Project vector against a standard vector with random structure.
Rand = Zero
do i=1,nA
  Rand = Rand+A(i)*Random_Molcas(iSeed)
end do
if (Rand < Zero) then
  A(:) = -A(:)
  Rand = -Rand
end if

end subroutine VecPhase
