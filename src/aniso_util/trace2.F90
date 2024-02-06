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

complex*16 function trace2(n,A,B)

use Constants, only: cZero

implicit none
! size of the square matrices A(n,n) and B(n,n)
integer, intent(in) :: n
complex(kind=8), intent(in) :: A(n,n), B(n,n)
! local variables
integer :: i, k

trace2 = cZero
do i=1,n
  do k=1,n
    trace2 = trace2+A(k,i)*B(i,k)
  end do
end do

end function trace2
