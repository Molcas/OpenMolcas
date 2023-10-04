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

subroutine ppgs_cvb(cvb)

use casvb_global, only: i2s_fr, ifnss1, mnion_fr, nconf_fr, nel_fr, nfrag, nS_fr, nvb_fr, vbdet
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
real(kind=wp) :: cvb(nvb)
integer(kind=iwp) :: icoffs_nconf, ifrag, ioffs_cvb, iS, kbasiscvb_kp, nelsing

! First applicable configuration with first possible spin in
! each fragment is set to perfect-pairing:
call dcopy_(nvb,[1.0e-2_wp],0,cvb,1)
ioffs_cvb = 0
icoffs_nconf = 0
do ifrag=1,nfrag
  nelsing = nel_fr(ifrag)-2*mnion_fr(ifrag)
  do iS=1,nS_fr(ifrag)
    if (i2s_fr(iS,ifrag) <= nelsing) then
      cvb(ifnss1(nelsing,i2s_fr(iS,ifrag))+ioffs_cvb) = One
      exit
    end if
  end do
  ioffs_cvb = ioffs_cvb+nvb_fr(ifrag)
  icoffs_nconf = icoffs_nconf+nconf_fr(ifrag)
end do
kbasiscvb_kp = kbasiscvb
kbasiscvb = 1
call str2vbc_cvb(cvb,vbdet)
kbasiscvb = kbasiscvb_kp
call vb2strc_cvb(vbdet,cvb)

return

end subroutine ppgs_cvb
