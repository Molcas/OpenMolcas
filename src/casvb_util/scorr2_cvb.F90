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

subroutine scorr2_cvb(cvbdet,dvbdet,evbdet,ssq,wvbdet,iperm)

use casvb_global, only: formAD, formAF
use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6

implicit none
#include "main_cvb.fh"
real(kind=wp) :: cvbdet(ndetvb), dvbdet(ndetvb), evbdet(ndetvb), ssq(norb,norb), wvbdet(ndetvb)
integer(kind=iwp) :: iperm(norb)
real(kind=wp), parameter :: cut = 1.0e-10_wp
integer(kind=iwp) :: i, mu, nu
real(kind=wp) :: phase, rsum, scheck, snorm, ssnorm, ssum, stot, tot
real(kind=wp), external :: ddot_

write(u6,'(/,1x,a)') 'Expectation values of (s(i)+s(j))**2'
snorm = ddot_(ndetvb,cvbdet,1,dvbdet,1)
ssnorm = ddot_(ndetvb,cvbdet,1,evbdet,1)
write(u6,formAF) ' Lower triangle uses SPIN function with Snorm=',ssnorm
write(u6,formAF) ' Upper triangle uses FULL function with Snorm=',snorm
! DLC
!snorm = one/snorm
!ssnorm = one/ssnorm
phase = (-one)**abs(nalf-nbet)
snorm = phase/snorm
ssnorm = phase/ssnorm
!! DLC
call fzero(ssq,norb*norb)
tot = Zero
stot = Zero
do mu=1,norb
  do nu=mu+1,norb
    ! Apply s_mu x s_nu to the wavefunction
    do i=1,norb
      iperm(i) = i
    end do
    iperm(mu) = nu
    iperm(nu) = mu
    call fmove_cvb(cvbdet,wvbdet,ndetvb)
    call permvb_cvb(wvbdet,iperm)
    rsum = one-ddot_(ndetvb,wvbdet,1,dvbdet,1)*snorm
    ssum = one-ddot_(ndetvb,wvbdet,1,evbdet,1)*ssnorm
    tot = tot+rsum
    stot = stot+ssum
    ssq(mu,nu) = rsum
    ssq(nu,mu) = ssum
  end do
end do
call mxprint_cvb(ssq,norb,norb,0)
tot = tot+r3by4*real(norb-2*norb*(norb-1)/2,kind=wp)
stot = stot+r3by4*real(norb-2*norb*(norb-1)/2,kind=wp)
scheck = Half*real(abs(nalf-nbet),kind=wp)*(Half*real(abs(nalf-nbet),kind=wp)+one)
if ((abs(tot-scheck) > cut) .or. (abs(stot-scheck) > cut)) write(u6,formAD) 'WARNING: spins ',stot,tot,scheck

return

end subroutine scorr2_cvb
