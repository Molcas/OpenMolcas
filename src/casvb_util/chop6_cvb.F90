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

subroutine chop6_cvb()

use casvb_global, only: icase6, mxdav, release

implicit real*8(a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
logical done

if (release(6)) call mfreer_cvb(lp(1))
release(6) = .true.
release(7) = .false.
lp(1) = mstackr_cvb(0)

call setcnt2_cvb(6,0)
if (icase6 == 1) then
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
else if (icase6 == 2) then
  ! Overlap-based Davidson optimization:
  iremain = mavailr_cvb()
  maxdav = min(mxiter,nvb,mxdav)

  memwrk = ndetvb+5*norb*norb+3*ihlf_cvb(norb+2*norb*norb)
  done = .false.
  do idav=maxdav,1,-1
    ! NEED is approx req. memory:
    need = 2*nvb*idav+2*nvb+idav+1000+memwrk
    if (need < iremain) then
      done = .true.
      exit
    end if
  end do
  if (.not. done) then
    idav = 0
    if (nvb == 0) then
      need = 1000+memwrk
      if (need < iremain) done = .true.
    end if
  end if
  if (.not. done) then
    write(6,*) ' Not enough memory for Davidson!',need,iremain
    call abend_cvb()
  end if
  maxdav = idav

else if (icase6 == 3) then
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

  done = .false.
  do idav=maxdav,1,-1
    ! NEED is approx req. memory:
    need = 3*nvb*idav+nvb+idav*(2*idav+3)+1000+mem_applyh+memwrk
    if (need < iremain) then
      done = .true.
      exit
    end if
  end do
  if (.not. done) then
    idav = 0
    if (nvb == 0) then
      need = 1000+memwrk
      if (need < iremain) done = .true.
    end if
  end if
  if (.not. done) then
    write(6,*) ' Not enough memory for Davidson!',need,iremain
    call abend_cvb()
  end if
  maxdav = idav

else if (icase6 == 4) then
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

end subroutine chop6_cvb
