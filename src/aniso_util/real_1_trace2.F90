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

function real_1_trace2(n,A)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: real_1_trace2
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: A(n,n)
integer(kind=iwp) :: i

real_1_trace2 = Zero
do i=1,n
  real_1_trace2 = real_1_trace2+A(i,i)
end do
real_1_trace2 = real_1_trace2/real(n,kind=wp)

end function real_1_trace2
