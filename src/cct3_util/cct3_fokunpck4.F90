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

subroutine cct3_fokunpck4(fok,fii,dimfok,dimfi)
! this routine distributes (Fok - dp) to Fii
! fok    - Fok matrix (I)
! fii    - Fii matrix (O)
! dimfok - dimension for Fok matrix - norb (I)
! dimfi  - dimension of occupied - no (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimfok, dimfi
real(kind=wp), intent(in) :: fok(dimfok,dimfok)
real(kind=wp), intent(out) :: fii(dimfi,dimfi)

!1 distribute Fok to Fii
fii(:,:) = fok(1:dimfi,1:dimfi)

return

end subroutine cct3_fokunpck4
