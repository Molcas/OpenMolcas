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
integer(kind=iwp) :: dimp, dimq
real(kind=wp) :: a(dimp,dimq**2), b(dimp,dimq,dimq)
integer(kind=iwp) :: p, q, qr, r
real(kind=wp) :: scalar

if (dimq > 1) then

  qr = 0
  do q=2,dimq
    do r=1,q-1
      qr = qr+1

      do p=1,dimp
        scalar = a(p,qr)
        b(p,q,r) = scalar
        b(p,r,q) = -scalar
      end do
    end do
  end do

end if

do q=1,dimq
  do p=1,dimp
    b(p,q,q) = Zero
  end do
end do

return

end subroutine expand3
