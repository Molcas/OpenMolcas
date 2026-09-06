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

!#define _RANDOM_CODE_
subroutine VecPhase(A,nA)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nA
real(kind=wp), intent(inout) :: A(nA)
integer(kind=iwp) :: i
real(kind=wp) :: Phase
#ifdef _RANDOM_CODE_
integer(kind=iwp) :: iSeed = 17
real(kind=wp), external :: Random_Molcas
#endif

Phase = Zero
do i=1,nA
# ifdef _RANDOM_CODE_
  ! Project vector against a standard vector with random structure.
  Phase = Phase+A(i)*Random_Molcas(iSeed)
# else
  Phase = Phase+A(i)*i
# endif
end do
if (Phase < Zero) A(:) = -A(:)

end subroutine VecPhase
