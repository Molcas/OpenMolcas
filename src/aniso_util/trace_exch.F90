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

complex*16 function trace_exch(n1,n2,A,B)

use Constants, only: cZero

implicit none
! size of the square matrices A(n1,n1,n2,n2) and B(n1,n1,n2,n2)
integer, intent(in) :: n1, n2
complex(kind=8), intent(in) :: A(n1,n1,n2,n2), B(n1,n1,n2,n2)
! local variables
integer :: i1, i2, k1, k2

trace_exch = cZero
do i1=1,n1
  do k1=1,n1
    do i2=1,n2
      do k2=1,n2
        trace_exch = trace_exch+A(i1,k1,i2,k2)*B(k1,i1,k2,i2)
      end do
    end do
  end do
end do

end function trace_exch
