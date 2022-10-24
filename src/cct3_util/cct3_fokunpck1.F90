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

subroutine cct3_fokunpck1(fok,dp,dimfok)
! this routine does Fok = Fok - dp
! fok    - Fok matrix (I/O)
! dp     - Diagonal part vector (I)
! dimfok - dimension for Fok matrix - norb (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimfok
real(kind=wp), intent(inout) :: fok(dimfok,dimfok)
real(kind=wp), intent(in) :: dp(dimfok)
integer(kind=iwp) :: p

!1 substract dp from Fok
do p=1,dimfok
  fok(p,p) = fok(p,p)-dp(p)
end do

return

end subroutine cct3_fokunpck1
