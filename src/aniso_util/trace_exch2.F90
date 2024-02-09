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

function trace_exch2(n1,n2,A,O1,O2)

use Constants, only: cZero
use Definitions, only: wp, iwp

implicit none
complex(kind=wp) :: trace_exch2
integer(kind=iwp), intent(in) :: n1, n2
complex(kind=wp), intent(in) :: A(n1,n1,n2,n2), O1(n1,n1), O2(n2,n2)
integer(kind=iwp) :: i1, i2, k1, k2

trace_exch2 = cZero
do i1=1,n1
  do k1=1,n1
    do i2=1,n2
      do k2=1,n2
        trace_exch2 = trace_exch2+A(i1,k1,i2,k2)*O1(k1,i1)*O2(k2,i2)
      end do
    end do
  end do
end do

end function trace_exch2
