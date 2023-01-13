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

subroutine cct3_fokunpck3(fok,fai,dimfok,dimfa,dimfi)
! this routine distributes (Fok - dp) to Fai
! fok    - Fok matrix (I)
! fai    - Fai matrix (O)
! dimfok - dimension for Fok matrix - norb (I)
! dimfa  - dimension of virtuals - nv (I)
! dimfi  - dimension of occupied - no (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimfok, dimfa, dimfi
real(kind=wp), intent(in) :: fok(dimfok,dimfok)
real(kind=wp), intent(out) :: fai(dimfa,dimfi)

!1 distribute Fok to Fai
fai(:,:) = fok(dimfi+1:dimfi+dimfa,1:dimfi)

return

end subroutine cct3_fokunpck3
