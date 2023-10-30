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

subroutine makecivbhs_cvb(civbh,civbs,orbs)
! Construct CIVBS ( = T(s) * CIVB ) & CIVBH ( = T(O)*H*T(O) * CIVB ):

use casvb_global, only: icnt_ci, ndet, norb
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: civbh(0:ndet), civbs(0:ndet)
real(kind=wp), intent(in) :: orbs(norb,norb)
integer(kind=iwp) :: icivbh, icivbs

icivbh = nint(civbh(0))
icivbs = nint(civbs(0))
if ((icnt_ci(icivbs) == 4) .and. (icnt_ci(icivbh) == 5)) then
  return
else if (icnt_ci(icivbs) == 4) then
  call applyth_cvb(civbh,orbs)
else if (icnt_ci(icivbs) == 5) then
  call applyts_cvb(civbs,orbs)
else
  call applyths_cvb(civbh,civbs,orbs)
end if
icnt_ci(icivbs) = 4
icnt_ci(icivbh) = 5

return

end subroutine makecivbhs_cvb
