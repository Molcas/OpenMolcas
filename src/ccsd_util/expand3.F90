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

subroutine expand3(a,b,dimp,dimq)
! expand a(p,qr) -> b(p,q,r)
! assumption: q>r, a(p,q,r)=-a(p,r,q)
! RISC version

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimp, dimq
real(kind=wp), intent(in) :: a(dimp,dimq**2)
real(kind=wp), intent(out) :: b(dimp,dimq,dimq)
integer(kind=iwp) :: q, qr, r

if (dimq > 1) then

  qr = 0
  do q=2,dimq
    do r=1,q-1
      qr = qr+1
      b(:,q,r) = a(:,qr)
      b(:,r,q) = -a(:,qr)
    end do
  end do

end if

do q=1,dimq
  b(:,q,q) = Zero
end do

return

end subroutine expand3
