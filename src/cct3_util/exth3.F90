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

subroutine exth3(a,b,dimp,dimq,dimr,q,nfact)
! this routine extracts A(p,q,r) -> B_q(p,r)
!
! a     - matrix a (Input)
! b     - matrix b (Output)
! dimp  - dimension of p (Input)
! dimq  - dimension of q (Input)
! dimr  - dimension of r (Input)
! q     - value of index q (Input)
! nfact - sign (+-1,0) (Input)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimp, dimq, dimr, q, nfact
real(kind=wp), intent(in) :: a(dimp,dimq,dimr)
real(kind=wp), intent(out) :: b(dimp,dimr)

if (nfact == 1) then
  b(:,:) = a(:,q,:)
else if (nfact == -1) then
  b(:,:) = -a(:,q,:)
else if (nfact == 0) then
  b(:,:) = Zero
end if

return

end subroutine exth3
