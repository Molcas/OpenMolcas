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

subroutine orthcvb_cvb(c,nparm1)

use casvb_global, only: cvbnrm_fr, nfrag, nvb_fr
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: c(*)
integer(kind=iwp) :: nparm1
#include "main_cvb.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: ifr_off, ifrag, ioffs
real(kind=wp), external :: ddot_

ioffs = nparm1-nprvb+1
if (nfrag <= 1) then
  call daxpy_(nprvb,-ddot_(nprvb,work(lv(2)),1,c(ioffs),1)/cvbnrm,work(lv(2)),1,c(ioffs),1)
else
  ifr_off = 0
  do ifrag=1,nfrag
    call daxpy_(nvb_fr(ifrag),-ddot_(nvb_fr(ifrag),work(ifr_off+lv(2)),1,c(ifr_off+ioffs),1)/cvbnrm_fr(ifrag),work(ifr_off+lv(2)), &
                1,c(ifr_off+ioffs),1)
    ifr_off = ifr_off+nvb_fr(ifrag)
  end do
end if

return

end subroutine orthcvb_cvb
