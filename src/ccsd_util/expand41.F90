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

subroutine expand41(a,b,dimpq,dimr,dims,dimp)
! expand a(pq,r,s) -> - b(p,q,s,r)
! assumption: p>q + antisymmetry
! RISC version

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimpq, dimr, dims, dimp
real(kind=wp), intent(in) :: a(dimpq,dimr,dims)
real(kind=wp), intent(out) :: b(dimp,dimp,dims,dimr)
integer(kind=iwp) :: p, pq, q, r, s
real(kind=wp) :: scalar

if (dimp > 1) then

  do s=1,dims
    do r=1,dimr

      pq = 0
      do p=2,dimp
        do q=1,p-1
          pq = pq+1

          scalar = a(pq,r,s)
          b(p,q,s,r) = -scalar
          b(q,p,s,r) = scalar

        end do
      end do
    end do
  end do

end if

do p=1,dimp
  b(p,p,:,:) = Zero
end do

return

end subroutine expand41
