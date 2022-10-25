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

subroutine cct3_expand40(a,b,dimpq,dimrs,dimp,dimr)
! expand a(pq,rs) -> b(p,q,r,s)
! assumption: p>q,r>s, + antisymmetry
! RISC version

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimpq, dimrs, dimp, dimr
real(kind=wp), intent(in) :: a(dimpq,dimrs)
real(kind=wp), intent(out) :: b(dimp,dimp,dimr,dimr)
integer(kind=iwp) :: p, pq, r, rs

if (dimp > 1) then

  rs = 0
  do r=1,dimr
    pq = 0
    do p=2,dimp
      b(p,1:p-1,r,1:r-1) = a(pq+1:pq+p-1,rs+1:rs+r-1)
      b(p,1:p-1,1:r-1,r) = -a(pq+1:pq+p-1,rs+1:rs+r-1)
      b(1:p-1,p,r,1:r-1) = -a(pq+1:pq+p-1,rs+1:rs+r-1)
      b(1:p-1,p,1:r-1,r) = a(pq+1:pq+p-1,rs+1:rs+r-1)
      pq = pq+p-1
    end do
    rs = rs+r-1
    b(:,:,r,r) = Zero
  end do

end if

do p=1,dimp
  b(p,p,:,:) = Zero
end do

return

end subroutine cct3_expand40
