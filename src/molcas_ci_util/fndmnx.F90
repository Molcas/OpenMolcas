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

real*8 function FNDMNX(VECTOR,NDIM,MINMAX)
! FIND SMALLEST(MINMAX=1) OR LARGEST(MINMAX=2)
! ABSOLUTE VALUE OF ELEMENTS IN VECTOR

implicit real*8(A-H,O-Z)
dimension VECTOR(*)

! jwk-cleanup
result = 0.0d0
if (MINMAX == 1) then
  result = abs(VECTOR(1))
  do I=2,NDIM
    result = min(result,abs(VECTOR(I)))
  end do
end if

if (MINMAX == 2) then
  result = abs(VECTOR(1))
  do I=2,NDIM
    result = max(result,abs(VECTOR(I)))
  end do
end if

FNDMNX = result

return

end function FNDMNX
