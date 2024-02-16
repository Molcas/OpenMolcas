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

subroutine dirvect(P1,R1,P2,R2,vec,dist)
! this Subroutine computes the directional vector of the origins of two points P1 and P2

use Constants, only: Zero

implicit none
real(kind=8), intent(in) :: P1(3) ! coords of the first point
real(kind=8), intent(in) :: P2(3) ! coords of the second point
real(kind=8), intent(in) :: R1(3,3) ! rot. matrix of the  first point to the general coordinate system
real(kind=8), intent(in) :: R2(3,3) ! rot. matrix of the second point to the general coordinate system
real(kind=8), intent(out) :: vec(3)
real(kind=8), intent(out) :: dist
! local variables
integer :: i, j
real(kind=8) :: C1(3), C2(3)
real(kind=8) :: distance
external :: distance

C1(:) = Zero
C2(:) = Zero

do i=1,3
  do j=1,3
    C1(i) = C1(i)+P1(j)*R1(i,j)
    C2(i) = C2(i)+P2(j)*R2(i,j)
  end do
end do
dist = distance(3,C1,C2)
do i=1,3
  vec(i) = (C1(i)-C2(i))/dist
end do

return

end subroutine dirvect
