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

subroutine fx_svb1_cvb(fx,fast,orbstry,cvbtry,civec,civecp,civbs,civb,cvbdet)

use casvb_global, only: formE, gjorb, ipr, memplenty, ndet, ndetvb, norb, nvb, ovraa_try, ovrab_try, svb
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(out) :: fx, cvbdet(ndetvb)
logical(kind=iwp), intent(in) :: fast
real(kind=wp), intent(in) :: orbstry(norb,norb)
real(kind=wp), intent(inout) :: cvbtry(nvb), civec(0:ndet), civecp(0:ndet), civbs(0:ndet), civb(0:ndet)

call str2vbc_cvb(cvbtry,cvbdet)
if (fast) then
  call makecivb_cvb(civec,civb,cvbdet,orbstry,cvbtry,1)
  call gaussj_cvb(orbstry,gjorb)
  call applyt_cvb(civb,gjorb)
  call proj_cvb(civb)
  call cinorm_cvb(civb,ovraa_try)
  if (memplenty) then
    call cidot_cvb(civec,civb,ovrab_try)
  else
    call cird_cvb(civecp,61001.2_wp)
    call cidot_cvb(civecp,civb,ovrab_try)
  end if
else
  call makecivb_cvb(civec,civb,cvbdet,orbstry,cvbtry,0)
  call makecivecp_cvb(civec,civecp,orbstry)
  call makecivbs_cvb(civbs,orbstry,cvbdet)

  call pvbdot_cvb(civb,civbs,ovraa_try)
  call pvbdot_cvb(civb,civecp,ovrab_try)
end if
svb = ovrab_try/sqrt(ovraa_try)
fx = svb
if (fast .and. (ipr(3) >= 2)) write(u6,formE) ' Svb :      ',svb

return

end subroutine fx_svb1_cvb
