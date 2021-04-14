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

integer function IDMIN(N,X,INCX)

! FINDS THE INDEX OF ELEMENT HAVING MIN. VALUE.

integer N, INCX
real*8 X(*)
integer I, IXX
real*8 MIN

IDMIN = 0
if (N >= 1) then
  IDMIN = 1
  if (N /= 1) then
    if (INCX /= 1) then
      ! CODE FOR INCREMENT NOT EQUAL TO 1
      IXX = 1
      MIN = X(1)
      IXX = IXX+INCX
      do I=2,N
        if (X(IXX) < MIN) then
          IDMIN = I
          MIN = X(IXX)
        end if
        IXX = IXX+INCX
      end do
    else
      ! CODE FOR INCREMENT EQUAL TO 1
      MIN = X(1)
      do I=2,N
        if (X(I) < MIN) then
          IDMIN = I
          MIN = X(I)
        end if
      end do
    end if
  end if
end if

return

end function IDMIN
