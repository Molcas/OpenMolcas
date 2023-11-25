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

use casvb_global, only: endvar, iciweights, icrit, ifhamil, ifinish, imethod, lcalccivbs, lcalcevb, lcalcsvb, lciweights, &
                        memplenty, nalf, nbet, nda, ndb, ndet, ndres, ndres_ok, norb, npcf, nv, variat
use Constants, only: Ten
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: icase, icritold, ifin_old, mavailr
real(kind=wp) :: strtcas
logical(kind=iwp) :: changed
logical(kind=iwp), external :: chpcmp_cvb, &  ! ... Change of dimensioning variables ...
                               ifcasci_cvb, & ! ... Files available ...
                               up2date_cvb    ! ... Make: up to date? ...

! This was in a common block, but never initialized or used elsewhere
strtcas = 0

changed = .false.
! CI vectors
call icomb_cvb(norb,nalf,nda)
call icomb_cvb(norb,nbet,ndb)
ndet = nda*ndb
ndres = 1+ndet

call mma_maxDBLE(mavailr)
memplenty = (mavailr > 9*ndet)

if (chpcmp_cvb(ndres)) changed = .true.
ndres_ok = (.not. changed) .and. (up2date_cvb('MEM3'))

if (changed) call touch_cvb('RDCAS')
if (chpcmp_cvb(nint(strtcas*Ten))) call touch_cvb('RDCAS')
call chpcmp2_cvb(icrit,icritold)
call chpcmp2_cvb(ifinish,ifin_old)
if (.not. ((icritold == 1) .and. (ifin_old == 0))) call touch_cvb('RDCAS')

if ((ifinish == 1) .or. (ifinish == 2)) then
  ! Logical variables for wavefunction analysis:
  ! Always evaluate svb/evb when possible
  !   (a) to get exact values (wfn is updated by optim in last iteration)
  !   (b) to generate ESYM (variational calculations)
  lcalcsvb = ifcasci_cvb()
  lcalcevb = ifhamil
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

end subroutine change4_cvb
