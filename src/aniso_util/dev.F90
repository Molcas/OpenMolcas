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

real*8 function dev(N,Fcal,Fexp)
! this function Returns the standard deviation between experimental
! data points and the computed ones
!     N --- number of data points ( Integer, input)
!  Fcal --- calculated array of size (N), Real(kind=8) ::, input
!  Fexp --- experimental array of size (N), Real(kind=8) ::, input
!   dev --- standard deviation, Real(kind=8) ::, output;

use Constants, only: Zero
use Definitions, only: wp

implicit none
integer, intent(in) :: N
real(kind=8), intent(in) :: Fcal(N), Fexp(N)
integer :: i
real(kind=8) :: diff, X

X = Zero
do i=1,N
  diff = Fcal(i)-Fexp(i)
  X = X+diff**2
end do
dev = sqrt(X/real(N,kind=wp))

return

end function dev
