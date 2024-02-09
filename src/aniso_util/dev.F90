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

function dev(N,Fcal,Fexp)
! this function Returns the standard deviation between experimental
! data points and the computed ones
!     N --- number of data points ( Integer, input)
!  Fcal --- calculated array of size (N), Real(kind=wp) ::, input
!  Fexp --- experimental array of size (N), Real(kind=wp) ::, input
!   dev --- standard deviation, Real(kind=wp) ::, output;

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: dev
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(in) :: Fcal(N), Fexp(N)
integer(kind=iwp) :: i
real(kind=wp) :: diff, X

X = Zero
do i=1,N
  diff = Fcal(i)-Fexp(i)
  X = X+diff**2
end do
dev = sqrt(X/real(N,kind=wp))

return

end function dev
