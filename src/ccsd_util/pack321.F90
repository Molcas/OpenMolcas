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

subroutine pack321(ap,am,b,dimp,dimq,dimr,rc)
! this routine does: B(p,q,r) = A+(p,q,r) - A-(p,r,q) for symq>symr

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimp, dimq, dimr
real(kind=wp), intent(in) :: ap(dimp,dimq,dimr), am(dimp,dimr,dimq)
real(kind=wp), intent(out) :: b(dimp,dimq,dimr)
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp) :: r

rc = 0
do r=1,dimr
  b(:,:,r) = ap(:,:,r)-am(:,r,:)
end do

return

end subroutine pack321
