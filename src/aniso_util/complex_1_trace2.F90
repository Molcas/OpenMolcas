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

complex*16 function complex_1_trace2(n,A)

use Constants, only: Zero, cOne

implicit none
! size of the square matrices A(n,n)
integer, intent(in) :: n
complex(kind=8), intent(in) :: A(n,n)
! local variables
integer :: i
complex(kind=8) :: FACT

complex_1_trace2 = Zero
FACT = n*cOne
do i=1,n
  complex_1_trace2 = complex_1_trace2+A(i,i)/FACT
end do

end function complex_1_trace2
