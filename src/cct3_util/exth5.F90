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

subroutine exth5(a,b,dimp,dimq,dimqr,q)
! this routine extracts A(p,qr) -> B_q(p,r)
!
! a     - matrix a (Input)
! b     - matrix b (Output)
! dimp  - dimension of p (Input)
! dimq  - dimension of q (and also r) (Input)
! dimqr - dimension of qr (Input)
! q     - value of index q (Input)

use CCT3_global, only: nshf
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimp, dimq, dimqr, q
real(kind=wp), intent(in) :: a(dimp,dimqr)
real(kind=wp), intent(out) :: b(dimp,dimq)
integer(kind=iwp) :: qr0, r, rq

if (q == 0) return

! r>q part
if (q > 1) then
  qr0 = nshf(q)
  b(:,1:q-1) = a(:,qr0+1:qr0+q-1)
end if

! r=q part
b(:,q) = Zero

! r<p part
if (q < dimq) then
  do r=q+1,dimq
    rq = nshf(r)+q
    b(:,r) = -a(:,rq)
  end do
end if

return

end subroutine exth5
