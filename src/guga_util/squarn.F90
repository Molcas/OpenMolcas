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

subroutine SQUARN(A,B,N)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: A(*)
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(out) :: B(N,N)
integer(kind=iwp) :: I, IIN

IIN = 2
do I=2,N
  B(I,1:I-1) = -A(IIN:IIN+I-2)
  B(1:I-1,I) = A(IIN:IIN+I-2)
  IIN = IIN+I
end do
call DCOPY_(N,[Zero],0,B,N+1)

return

end subroutine SQUARN
