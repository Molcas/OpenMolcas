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

real*8 function real_1_trace2(n,A)

implicit none
integer, parameter :: wp = kind(0.d0)
! size of the square matrices A(n,n)
integer, intent(in) :: n
real(kind=8), intent(in) :: A(n,n)
! local variables
integer :: i

real_1_trace2 = 0.0_wp
do i=1,n
  real_1_trace2 = real_1_trace2+A(i,i)/dble(n)
end do

end function real_1_trace2
