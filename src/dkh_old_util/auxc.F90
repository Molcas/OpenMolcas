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

subroutine AUXC(N,NBC,TC,ANSC)

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N, NBC
real(kind=wp), intent(in) :: TC
real(kind=wp), intent(out) :: ANSC
integer(kind=iwp) :: I
real(kind=wp) :: DUM, X, Y

X = One/(One+TC)
Y = TC*X
X = sqrt(X**(NBC+1))
DUM = One
if (N > 1) then
  do I=N,2,-1
    DUM = DUM*Y*real(NBC+2*I-3,kind=wp)/real(2*I-2,kind=wp)+One
  end do
end if
ANSC = X*DUM

return

end subroutine AUXC
