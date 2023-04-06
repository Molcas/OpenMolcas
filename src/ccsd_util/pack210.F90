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

subroutine pack210(a,b,dimpq,dimp,rc)
! this routine does: B(pq) = A(p,q) - A(q,p) for symp=symq

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimpq, dimp
real(kind=wp), intent(in) :: a(dimp,dimp)
real(kind=wp), intent(inout) :: b(dimpq)
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp) :: p, pq, q

rc = 0
if (dimp > 1) then

  pq = 0
  do p=2,dimp
    do q=1,p-1
      pq = pq+1
      b(pq) = a(p,q)-a(q,p)
    end do
  end do

else
  ! RC=1 : dimp is less than 2
  rc = 1
end if

return

end subroutine pack210
