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

use casvb_global, only: cvb, cvbnrm, cvbnrm_fr, nfrag, nprvb, nvb_fr
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: c(*)
integer(kind=iwp), intent(in) :: nparm1
integer(kind=iwp) :: ifr_off, ifrag, ioff2, ioffs, nprvb2
real(kind=wp) :: f
real(kind=wp), external :: ddot_

ioffs = nparm1-nprvb
if (nfrag <= 1) then
  f = ddot_(nprvb,cvb(1:nprvb),1,c(ioffs+1:ioffs+nprvb),1)/cvbnrm
  c(ioffs+1:ioffs+nprvb) = c(ioffs+1:ioffs+nprvb)-f*cvb(1:nprvb)
else
  ifr_off = 0
  do ifrag=1,nfrag
    nprvb2 = nvb_fr(ifrag)
    ioff2 = ifr_off+ioffs
    f = ddot_(nprvb2,cvb(ifr_off+1:ifr_off+nprvb2),1,c(ioff2+1:ioff2+nprvb2),1)/cvbnrm_fr(ifrag)
    c(ioff2+1:ioff2+nprvb2) = c(ioff2+1:ioff2+nprvb2)-f*cvb(ifr_off+1:ifr_off+nprvb2)
    ifr_off = ifr_off+nprvb2
  end do
end if

return

end subroutine orthcvb_cvb
