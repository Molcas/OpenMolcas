!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine applyts_cvb(civbs,orbs)
! Apply T(s) to CIVBS:

use casvb_global, only: gjorb, gjorb2, gjorb3, ndet, norb, proj
use Definitions, only: wp

implicit none
real(kind=wp), intent(inout) :: civbs(0:ndet)
real(kind=wp), intent(in) :: orbs(norb,norb)

call makegjorbs_cvb(orbs)

if (.not. proj) then
  call applyt_cvb(civbs,gjorb3)
else
  call applyt_cvb(civbs,gjorb)
  call proj_cvb(civbs)
  call applyt_cvb(civbs,gjorb2)
end if

return

end subroutine applyts_cvb
