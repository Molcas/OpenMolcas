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
integer(kind=iwp) :: q, qr

qr = 0
do q=1,dimq
  b(:,q,1:q-1,:) = a(:,qr+1:qr+q-1,:)
  b(:,1:q-1,q,:) = -a(:,qr+1:qr+q-1,:)
  qr = qr+q-1
  b(:,q,q,:) = Zero
end do

return

end subroutine cct3_expand2
