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

subroutine expand4_12(AA,BB,d1,d2,d3)
! this routine does:
!
! A(ab,i,j) -> A(a,b,i,j)

implicit none
integer d1, d2, d3, a, b, i, j, ab
real*8 AA(1:(d1*(d1+1))/2,d2,d3), BB(1:d1,1:d1,1:d2,1:d3)

ab = 0
do a=1,d1
  do b=1,a
    ab = ab+1
    do i=1,d2
      do j=1,d3
        BB(a,b,i,j) = AA(ab,i,j)
        if (a /= b) BB(b,a,j,i) = AA(ab,i,j)
      end do
    end do
  end do
end do

return

end subroutine expand4_12
