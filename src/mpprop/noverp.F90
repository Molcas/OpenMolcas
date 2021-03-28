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

subroutine NoverP(n,i,x)

implicit real*8(a-h,o-z)

rn = 1.0d0
rp = 1.0d0
if ((i == 0) .or. (i == n)) then
  x = 1
else
  do j=1,i
    rn = rn*(n-j+1)
    rp = rp*j
  end do
  x = rn/rp
end if

return

end subroutine NoverP
