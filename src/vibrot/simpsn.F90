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
!
! THIS SUBROUTINE EVALUATES THE INTEGRAL INT(XDX) USING SIMPSON'S
! RULE. THE NUMBER OF STEPS IS NDIM (ASSUMED ODD) WITH THE GRID
! SIZE H. THE RESULTING VALUE OF THE INTEGRAL IS RETURNED IN F.
!
! ********** RELEASE 81 04 09 **********

subroutine SIMPSN(X,H,NDIM,F)

use Constants, only: Two, Six
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: X(*), H
integer(kind=iwp), intent(in) :: NDIM
real(kind=wp), intent(out) :: F
integer(kind=iwp) :: I, NDIM1

F = X(1)+X(NDIM)
NDIM1 = NDIM-1
do I=2,NDIM1
  F = F+Two*X(I)
  if (mod(I,2) == 0) F = F+Two*X(I)
end do
F = F*H/Six

return

end subroutine SIMPSN
