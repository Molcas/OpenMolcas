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

subroutine expand40(a,b,dimpq,dimrs,dimp,dimr)
! expand a(pq,rs) -> b(p,q,r,s)
! assumption: p>q,r>s, + antisymmetry
! RISC version

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimpq, dimrs, dimp, dimr
real(kind=wp), intent(in) :: a(dimpq,dimrs)
real(kind=wp), intent(out) :: b(dimp,dimp,dimr,dimr)
integer(kind=iwp) :: p, pq, q, r, rs, s
real(kind=wp) :: scalar

if ((dimp > 1) .and. (dimr > 1)) then

  rs = 0
  do r=2,dimr
    do s=1,r-1
      rs = rs+1

      pq = 0
      do p=2,dimp
        do q=1,p-1
          pq = pq+1

          scalar = a(pq,rs)
          b(p,q,r,s) = scalar
          b(p,q,s,r) = -scalar
          b(q,p,r,s) = -scalar
          b(q,p,s,r) = scalar

        end do
      end do
    end do
  end do

end if

do r=1,dimr
  b(:,:,r,r) = Zero
end do

do p=1,dimp
  b(p,p,:,:) = Zero
end do

return

end subroutine expand40
