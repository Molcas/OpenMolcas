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

subroutine fokunpck2(fok,faa,dimfok,dimfa,shift)
! this routine distributes (Fok - dp) to Faa
! fok    - Fok matrix (I)
! faa    - Faa matrix (O)
! dimfok - dimension for Fok matrix - norb (I)
! dimfa  - dimension of virtuals - nv (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: dimfok, dimfa, shift
real(kind=wp) :: fok(dimfok,dimfok), faa(dimfa,dimfa)

!1 distribute Fok to Faa
faa(1:dimfa,1:dimfa) = fok(shift+1:shift+dimfa,shift+1:shift+dimfa)

return

end subroutine fokunpck2
