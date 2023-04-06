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

subroutine pack320(a,b,dimp,dimqr,dimq,rc)
! this routine does: B(p,qr) = A(p,q,r) - A(p,r,q) for symq=symr

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimp, dimqr, dimq
real(kind=wp), intent(in) :: a(dimp,dimq,dimq)
real(kind=wp), intent(inout) :: b(dimp,dimqr)
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp) :: p, q, qr, r

rc = 0
if (dimq > 1) then

  qr = 0
  do q=2,dimq
    do r=1,q-1
      qr = qr+1
      do p=1,dimp
        b(p,qr) = a(p,q,r)-a(p,r,q)
      end do
    end do
  end do

else
  ! RC=1 : dimp is less than 2
  rc = 1
end if

return

end subroutine pack320
