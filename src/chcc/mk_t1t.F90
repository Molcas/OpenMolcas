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

subroutine Mk_T1t(T1,H,dimbepp,no,nv,addbepp)
! this routine does:
! H(i,be") <- T1o(be,i)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimbepp, no, nv, addbepp
real(kind=wp), intent(in) :: T1(nv,no)
real(kind=wp), intent(out) :: H(no,dimbepp)
integer(kind=iwp) :: bepp

do bepp=1,dimbepp
  H(:,bepp) = T1(addbepp+bepp,:)
end do

return

end subroutine Mk_T1t
