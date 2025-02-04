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

subroutine Gen_iSD4(iS,jS,kS,lS,iSD,nSD,iSD4)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: iS, jS, kS, lS, nSD, iSD(0:nSD,1024)
integer(kind=iwp), intent(out) :: iSD4(0:nSD,4)

iSD4(:,1) = iSD(:,iS)
iSD4(:,2) = iSD(:,jS)
iSD4(:,3) = iSD(:,kS)
iSD4(:,4) = iSD(:,lS)

end subroutine Gen_iSD4
