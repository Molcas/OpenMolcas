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

complex*16 function trace(n,A,B)

implicit none
integer, parameter :: wp = kind(0.d0)
! size of the square matrices A(n,n) and B(n,n)
integer, intent(in) :: n
complex(kind=8), intent(in) :: A(n,n), B(n,n)
! local variables
integer :: i, k

trace = (0.0_wp,0.0_wp)
do i=1,n
  do k=1,n
    trace = trace+A(i,k)*B(k,i)
  end do
end do

end function trace
