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

subroutine ISWAP(N,X,incX,Y,incY)
!***********************************************************************
!                                                                      *
!     Interchange vectors X and Y                                      *
!                                                                      *
!     calling arguments:                                               *
!     N       : Integer, input.                                        *
!               Number of input elements.                              *
!     X       : Array of Integer                                       *
!               Vector X                                               *
!     incX    : Integer, input.                                        *
!               Stride of vector X.                                    *
!     Y       : Array of Integer                                       *
!               Vector Y                                               *
!     incY    : Integer, input.                                        *
!               Stride of vector Y.                                    *
!                                                                      *
!***********************************************************************

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: N, incX, incY
integer(kind=iwp), intent(inout) :: X(*), Y(*)
integer(kind=iwp) :: i, iX, iY, Temp

if (N < 0) then
  write(u6,*)
  write(u6,*) '  *** Error in subroutine ISWAP ***'
  write(u6,*) '  Invalid number of elements in vectors X and Y :'
  write(u6,*) '  N must be larger than zero'
  write(u6,*)
  call Abend()
end if

iX = 1
if (incX < 0) iX = 1+(1-N)*incX
iY = 1
if (incY < 0) iY = 1+(1-N)*incY

do i=0,N-1
  Temp = Y(iY+i*incY)
  Y(iY+i*incY) = X(iX+i*incX)
  X(iX+i*incX) = Temp
end do

return

end subroutine ISWAP
