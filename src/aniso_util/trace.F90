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

function trace(n,A,B)

use Constants, only: cZero
use Definitions, only: wp, iwp

implicit none
complex(kind=wp) :: trace
integer(kind=iwp), intent(in) :: n
complex(kind=wp), intent(in) :: A(n,n), B(n,n)
integer(kind=iwp) :: i, k

trace = cZero
do i=1,n
  do k=1,n
    trace = trace+A(i,k)*B(k,i)
  end do
end do

end function trace
