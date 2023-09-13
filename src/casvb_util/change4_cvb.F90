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

subroutine change4_cvb()

implicit real*8(a-h,o-z)
logical changed
logical ndres_ok
! ... Files/Hamiltonian available ...
logical, external :: ifcasci_cvb, ifhamil_cvb
! ... Make: up to date? ...
logical, external :: up2date_cvb
! ... Change of dimensioning variables ...
logical, external :: chpcmp_cvb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
#include "casinfo_cvb.fh"
#include "rls_cvb.fh"
#include "WrkSpc.fh"
save ndres_ok

changed = .false.
! CI vectors
call icomb_cvb(norb,nalf,nda)
call icomb_cvb(norb,nbet,ndb)
ndet = nda*ndb
ndres = 3+ndet

memplenty = (mavailr_cvb() > 9*ndet)

if (chpcmp_cvb(ndres)) changed = .true.
ndres_ok = (.not. changed) .and. (up2date_cvb('MEM3'))

if (changed) call touch_cvb('RDCAS')
if (chpcmp_cvb(nint(strtcas*1d1))) call touch_cvb('RDCAS')
call chpcmp2_cvb(icrit,icritold)
call chpcmp2_cvb(ifinish,ifin_old)
if (.not. ((icritold == 1) .and. (ifin_old == 0))) call touch_cvb('RDCAS')

if ((ifinish == 1) .or. (ifinish == 2)) then
  ! Logical variables for wavefunction analysis:
  ! Always evaluate svb/evb when possible
  !   (a) to get exact values (wfn is updated by optim in last iteration)
  !   (b) to generate ESYM (variational calculations)
  lcalcsvb = ifcasci_cvb()
  lcalcevb = ifhamil_cvb()
  lcalccivbs = .true.
  lciweights = ((npcf > -2) .and. ((.not. variat) .or. endvar) .and. (iciweights > 0))
end if
if ((imethod == 11) .and. (ifinish == 0)) then
  nv = 2
  icase = 2
else if ((imethod /= 4) .and. (ifinish == 0)) then
  if (memplenty) then
    nv = 8
    icase = 1
  else
    if ((icrit == 2) .and. (imethod /= 6)) then
      ! No need for CITMP:
      nv = 3
      icase = 5
    else
      nv = 5
      icase = 2
    end if
  end if
else if ((imethod == 4) .and. (ifinish == 0)) then
  nv = 2
  icase = 2
else if ((imethod == 12) .and. (ifinish == 0)) then
  nv = 8
  icase = 7
else if ((ifinish == 1) .or. (ifinish == 2)) then
  nv = 5
  if (.not. lciweights) then
    nv = nv-1
    if ((.not. lcalcevb) .or. lcalccivbs) nv = nv-1
    if ((.not. lcalcsvb) .or. lcalccivbs) nv = nv-1
  end if
  icase = 3
else
  nv = 3
  icase = 4
end if
if (chpcmp_cvb(nv)) changed = .true.
if (chpcmp_cvb(icase)) changed = .true.
if (changed) call touch_cvb('MEM4')

return

entry chop4_cvb()
if (release(4)) call mfreer_cvb(lc(1)-2)
release(4) = .true.
release(5) = .false.

! CIVECP and CIVBH share memory --> LC(2)
do iv=1,nv
  lc(iv) = mstackr_cvb(ndres)
end do
do iv=1,nv
  work(lc(iv)) = zero
  work(lc(iv)+1) = zero
end do
do iv=1,nv
  lc(iv) = lc(iv)+2
end do
! Fix to put in "objects":
do iv=1,nv
  call creatci_cvb(iv,work(lc(iv)),lc(iv)+1,nint(work(lc(iv)-2)),work(lc(iv)-1))
  if (.not. ndres_ok) call setcnt_cvb(work(lc(iv)),0)
end do
!-- ins
if (((ifinish == 1) .or. (ifinish == 2)) .and. (.not. lciweights)) then
  if (((.not. lcalcevb) .or. lcalccivbs) .and. ((.not. lcalcsvb) .or. lcalccivbs)) then
    lc(3) = lc(1)
    lc(4) = lc(1)
  else if ((.not. lcalcevb) .or. lcalccivbs) then
    lc(4) = lc(3)
  else if ((.not. lcalcsvb) .or. lcalccivbs) then
    lc(4) = lc(3)
    lc(3) = lc(1)
  end if
end if
if (.not. memplenty) lc(nv+1) = lc(1)
if ((nv == 3) .or. (nv == 4) .or. (nv == 5)) then
  lc(6) = lc(2)
  lc(7) = lc(3)
  lc(8) = lc(4)
end if

return

end subroutine change4_cvb
