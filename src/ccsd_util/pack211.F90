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

subroutine pack211(ap,am,b,dimp,dimq,rc)
! this routine does: B(p,q) = A+(p,q) - A-(q,p) for symp>symq

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimp, dimq
real(kind=wp), intent(in) :: ap(dimp,dimq), am(dimq,dimp)
real(kind=wp), intent(out) :: b(dimp,dimq)
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp) :: q

rc = 0
do q=1,dimq
  b(:,q) = ap(:,q)-am(q,:)
end do

return

end subroutine pack211
