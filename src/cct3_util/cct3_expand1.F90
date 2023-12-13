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

subroutine cct3_expand1(a,b,dimpq,dimr,dimp)
! expand a(pq,r) -> b(p,q,r)
! assumption: p>q, a(p,q,r)=-a(q,p,r)
! RISC version

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimpq, dimr, dimp
real(kind=wp), intent(in) :: a(dimpq,dimr)
real(kind=wp), intent(out) :: b(dimp,dimp,dimr)
integer(kind=iwp) :: p, pq

pq = 0
do p=1,dimp
  b(p,1:p-1,:) = a(pq+1:pq+p-1,:)
  b(1:p-1,p,:) = -a(pq+1:pq+p-1,:)
  pq = pq+p-1
  b(p,p,:) = Zero
end do

return

end subroutine cct3_expand1
