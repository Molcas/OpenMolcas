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

function complex_1_trace2(n,A)

use Constants, only: Zero, cOne
use Definitions, only: wp, iwp

implicit none
complex(kind=wp) :: complex_1_trace2
integer(kind=iwp), intent(in) :: n
complex(kind=wp), intent(in) :: A(n,n)
integer(kind=iwp) :: i

complex_1_trace2 = Zero
do i=1,n
  complex_1_trace2 = complex_1_trace2+A(i,i)
end do
complex_1_trace2 = complex_1_trace2/(real(n,kind=wp)*cOne)

end function complex_1_trace2
