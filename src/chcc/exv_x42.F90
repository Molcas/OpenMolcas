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

subroutine ExV_X42(Vp,V,dimab,no)
! this routine does:
! Vp(a,b,i) <- V(a,b,i,i)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimab, no
real(kind=wp), intent(out) :: Vp(dimab,no)
real(kind=wp), intent(in) :: V(dimab,no,no)
integer(kind=iwp) :: i

do i=1,no
  Vp(:,i) = V(:,i,i)
end do

return

end subroutine ExV_X42
