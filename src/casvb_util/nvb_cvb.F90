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

function nvb_cvb(kbasis_loc)

use casvb_global, only: absym, ndetvb, ndetvb_fr, ndetvb2_fr, nfrag, nvb_fr, nvbr_fr
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: nvb_cvb
integer(kind=iwp), intent(in) :: kbasis_loc
integer(kind=iwp) :: ifrag, ndetvb2, nvb_loc, nvbr

ndetvb = 0
ndetvb2 = 0
nvbr = 0
do ifrag=1,nfrag
  if (kbasis_loc /= 6) then
    nvb_fr(ifrag) = nvbr_fr(ifrag)
  else
    if (absym(1)) then
      nvb_fr(ifrag) = ndetvb2_fr(ifrag)
    else
      nvb_fr(ifrag) = ndetvb_fr(ifrag)
    end if
  end if
  ndetvb = ndetvb+ndetvb_fr(ifrag)
  ndetvb2 = ndetvb2+ndetvb2_fr(ifrag)
  nvbr = nvbr+nvbr_fr(ifrag)
end do

if (kbasis_loc /= 6) then
  nvb_loc = nvbr
else
  if (absym(1)) then
    nvb_loc = ndetvb2
  else
    nvb_loc = ndetvb
  end if
end if
nvb_cvb = nvb_loc

return

end function nvb_cvb
