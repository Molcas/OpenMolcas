!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Anders Bernhardsson                                    *
!***********************************************************************

function Fact(R)

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: Fact
real(kind=wp), intent(in) :: R
integer(kind=iwp) :: i, j, n

n = nint(R)
i = 1
if (n == 0) then
  Fact = One
  return
end if
do j=1,n
  i = i*j
end do
Fact = real(i,kind=wp)

end function Fact
