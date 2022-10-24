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

subroutine cct3_expand2(a,b,dimp,dimqr,dims,dimq)
! expand a(p,qr,s) -> b(p,q,r,s)
! assumption: p>q, a(p,q,r,s)=-a(p,r,q,s)
! RISC version

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimp, dimqr, dims, dimq
real(kind=wp), intent(in) :: a(dimp,dimqr,dims)
real(kind=wp), intent(out) :: b(dimp,dimq,dimq,dims)
integer(kind=iwp) :: p, q, qr, r, s
real(kind=wp) :: scalar

if (dimq > 1) then

  do s=1,dims
    qr = 0
    do q=2,dimq
      do r=1,q-1
        qr = qr+1
        do p=1,dimp
          scalar = a(p,qr,s)
          b(p,q,r,s) = scalar
          b(p,r,q,s) = -scalar
        end do
      end do
    end do
  end do

end if

do s=1,dims
  do q=1,dimq
    do p=1,dimp
      b(p,q,q,s) = Zero
    end do
  end do
end do

return

end subroutine cct3_expand2
