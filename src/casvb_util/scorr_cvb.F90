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

subroutine scorr_cvb(cvbdet,dvbdet,evbdet)

use Index_Functions, only: nTri_Elem
use casvb_global, only: formAD, formAF, nalf, nbet, ndetvb, norb
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: cvbdet(ndetvb), dvbdet(ndetvb), evbdet(ndetvb)
real(kind=wp), parameter :: cut = 1.0e-10_wp
integer(kind=iwp) :: i, mu, nu
real(kind=wp) :: phase, rsum, scheck, snorm, ssnorm, ssum, stot, tot
integer(kind=iwp), allocatable :: iperm(:)
real(kind=wp), allocatable :: ssq(:,:), wvbdet(:)
real(kind=wp), external :: ddot_

call mma_allocate(ssq,norb,norb,label='ssq')
call mma_allocate(wvbdet,ndetvb,label='wvbdet')
call mma_allocate(iperm,norb,label='iperm')

write(u6,'(/,1x,a)') 'Expectation values of (s(i)+s(j))**2'
snorm = ddot_(ndetvb,cvbdet,1,dvbdet,1)
ssnorm = ddot_(ndetvb,cvbdet,1,evbdet,1)
write(u6,formAF) ' Lower triangle uses SPIN function with Snorm=',ssnorm
write(u6,formAF) ' Upper triangle uses FULL function with Snorm=',snorm
! DLC
!snorm = One/snorm
!ssnorm = One/ssnorm
phase = (-One)**abs(nalf-nbet)
snorm = phase/snorm
ssnorm = phase/ssnorm
!! DLC
ssq(:,:) = Zero
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
    wvbdet(:) = cvbdet(:)
    call permvb_cvb(wvbdet,iperm)
    rsum = One-ddot_(ndetvb,wvbdet,1,dvbdet,1)*snorm
    ssum = One-ddot_(ndetvb,wvbdet,1,evbdet,1)*ssnorm
    tot = tot+rsum
    stot = stot+ssum
    ssq(mu,nu) = rsum
    ssq(nu,mu) = ssum
  end do
end do
call mxprint_cvb(ssq,norb,norb,0)
tot = tot+0.75_wp*real(norb-2*nTri_Elem(norb-1),kind=wp)
stot = stot+0.75_wp*real(norb-2*nTri_Elem(norb-1),kind=wp)
scheck = Half*real(abs(nalf-nbet),kind=wp)*(Half*real(abs(nalf-nbet),kind=wp)+One)
if ((abs(tot-scheck) > cut) .or. (abs(stot-scheck) > cut)) write(u6,formAD) 'WARNING: spins ',stot,tot,scheck

call mma_deallocate(ssq)
call mma_deallocate(wvbdet)
call mma_deallocate(iperm)

return

end subroutine scorr_cvb
