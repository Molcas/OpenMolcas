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

subroutine exMap3_231(A,B,d1,d2)
! this routine does:
!
! A (a,bc) -> B(b,c,a)

implicit none
integer d1, d2, i1, i2, i3, i23
real*8 A(1:d1,1:(d2*(d2+1))/2)
real*8 B(1:d2,1:d2,1:d1)

i23 = 0
do i2=1,d2
  do i3=1,i2
    i23 = i23+1
    do i1=1,d1

      B(i2,i3,i1) = A(i1,i23)
      B(i3,i2,i1) = A(i1,i23)

    end do
  end do
end do

return

end subroutine exMap3_231
