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

subroutine pack310(a,b,dimpq,dimr,dimp,rc)
! this routine does: B(pq,r) = A(p,q,r) - A(q,p,r) for symp=symq

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimpq, dimr, dimp
real(kind=wp), intent(in) :: a(dimp,dimp,dimr)
real(kind=wp), intent(inout) :: b(dimpq,dimr)
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp) :: p, pq, q, r

rc = 0
if (dimp > 1) then

  do r=1,dimr
    pq = 0
    do p=2,dimp
      do q=1,p-1
        pq = pq+1
        b(pq,r) = a(p,q,r)-a(q,p,r)
      end do
    end do
  end do

else
  ! RC=1 : dimp is less than 2
  rc = 1
end if

return

end subroutine pack310
