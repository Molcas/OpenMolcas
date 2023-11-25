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

subroutine orthcvb_init_cvb()

use casvb_global, only: cvb, cvbnrm, cvbnrm_fr, nfrag, nvb, nvb_fr
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: ifr_off, ifrag
real(kind=wp), external :: ddot_

if (nfrag <= 1) then
  cvbnrm = ddot_(nvb,cvb,1,cvb,1)
else
  ifr_off = 1
  do ifrag=1,nfrag
    cvbnrm_fr(ifrag) = ddot_(nvb_fr(ifrag),cvb(ifr_off:),1,cvb(ifr_off:),1)
    ifr_off = ifr_off+nvb_fr(ifrag)
  end do
end if

return

end subroutine orthcvb_init_cvb
