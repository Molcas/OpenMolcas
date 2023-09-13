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

subroutine change6_cvb()

implicit real*8(a-h,o-z)
logical changed
! ... Change of dimensioning variables ...
logical, external :: chpcmp_cvb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
#include "rls_cvb.fh"
#include "davtune_cvb.fh"
save icase

changed = .false.
if (changed) call touch_cvb('CIFREE')

nprorb = norb*(norb-1)
if (strucopt) then
  nprvb = nvb
  npr = nprorb+nvb
else if (.not. ((imethod == 4) .or. (imethod == 6))) then
  nprvb = 0
  npr = nprorb
else
  npr = 0
end if
if (chpcmp_cvb(npr)) changed = .true.
if ((.not. ((imethod == 4) .or. (imethod == 6))) .and. (ifinish == 0)) then
  ! Standard 2nd-order procedure:
  icase = 1
else if ((imethod == 4) .and. (icrit == 1) .and. (ifinish == 0)) then
  ! Overlap-based Davidson
  icase = 2
else if ((imethod == 4) .and. (icrit == 2) .and. (ifinish == 0)) then
  ! Energy-based Davidson
  icase = 3
else if ((imethod == 6) .or. (ifinish == 1) .or. (ifinish == 2)) then
  ! No arrays needed
  icase = 4
else
  icase = 5
end if
if (chpcmp_cvb(icase)) changed = .true.
if (changed) call touch_cvb('MEM6')

return

entry chop6_cvb()
if (release(6)) call mfreer_cvb(lp(1))
release(6) = .true.
release(7) = .false.
lp(1) = mstackr_cvb(0)

call setcnt2_cvb(6,0)
if (icase == 1) then
  ! Standard non-linear optimization procedure:
  lp(1) = mstackr_cvb(norb*norb+nvb+1+mxirrep)
  lp(2) = mstackr_cvb(npr)
  lq(3) = mstackr_cvb(nprorb*nprorb)
  lq(4) = mstackr_cvb(norb**4)
  lp(5) = mstackr_cvb(npr)
  lp(6) = mstackr_cvb(npr)
  lq(7) = mstackr_cvb(npr)
  lq(8) = mstackr_cvb(npr)
  lq(9) = mstackr_cvb(norb*norb)
  ! Vec1 work array
  lq(10) = mstackr_cvb(max(npr,ndetvb))
else if (icase == 2) then
  ! Overlap-based Davidson optimization:
  iremain = mavailr_cvb()
  maxdav = min(mxiter,nvb,mxdav)

  memwrk = ndetvb+5*norb*norb+3*ihlf_cvb(norb+2*norb*norb)
  do idav=maxdav,1,-1
    ! NEED is approx req. memory:
    need = 2*nvb*idav+2*nvb+idav+1000+memwrk
    if (need < iremain) goto 2
  end do
  idav = 0
  if (nvb == 0) then
    need = 1000+memwrk
    if (need < iremain) goto 2
  end if
  write(6,*) ' Not enough memory for Davidson!',need,iremain
  call abend_cvb()
2 maxdav = idav

else if (icase == 3) then
  ! Energy-based Davidson optimization:
  iremain = mavailr_cvb()
  maxdav = min(mxiter,nvb,mxdav)

  mem_applyh = ndet+neread
  ncimx = 0
  do ir=1,nirrep
    ncimx = max(ncimx,ncivb(ir))
  end do
  if (ncimx /= ndet) mem_applyh = mem_applyh+ncimx
  memwrk = ndetvb+3*norb*norb+2*ihlf_cvb(norb+2*norb*norb)

  do idav=maxdav,1,-1
    ! NEED is approx req. memory:
    need = 3*nvb*idav+nvb+idav*(2*idav+3)+1000+mem_applyh+memwrk
    if (need < iremain) goto 12
  end do
  idav = 0
  if (nvb == 0) then
    need = 1000+memwrk
    if (need < iremain) goto 12
  end if
  write(6,*) ' Not enough memory for Davidson!',need,iremain
  call abend_cvb()
12 maxdav = idav

else if (icase == 4) then
  ! Wavefunction analysis:
  mstackr_cvb0 = mstackr_cvb(0)
  if (((.not. variat) .or. endvar) .and. ((ivbweights > 1) .or. (ishstruc == 1))) then
    lp(1) = mstackr_cvb(nvb*nvb)
    lp(2) = mstackr_cvb(nvb*nvb)
  else
    lp(1) = mstackr_cvb0
    lp(2) = mstackr_cvb0
  end if
  do i=3,11
    lp(i) = mstackr_cvb0
  end do
end if

return

end subroutine change6_cvb
