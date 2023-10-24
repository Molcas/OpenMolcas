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

subroutine chop7_cvb()

use casvb_global, only: cvbdet, cvbsspn, cvbstot, cvbtry, dvbdet, evbdet, gjorb, gjorb2, gjorb3, icase7, ndetvb, norb, nvb, &
                        orbinv, orbstry, owrk2, release, sorbs, vbdet
use stdalloc, only: mma_allocate, mma_deallocate

implicit none

if (release(7)) then
  if (allocated(orbinv)) call mma_deallocate(orbinv)
  if (allocated(sorbs)) call mma_deallocate(sorbs)
  if (allocated(owrk2)) call mma_deallocate(owrk2)
  if (allocated(gjorb%r)) call mma_deallocate(gjorb%r)
  if (allocated(gjorb%i1)) call mma_deallocate(gjorb%i1)
  if (allocated(gjorb%i2)) call mma_deallocate(gjorb%i2)
  if (allocated(gjorb2%r)) call mma_deallocate(gjorb2%r)
  if (allocated(gjorb2%i1)) call mma_deallocate(gjorb2%i1)
  if (allocated(gjorb2%i2)) call mma_deallocate(gjorb2%i2)
  if (allocated(gjorb3%r)) call mma_deallocate(gjorb3%r)
  if (allocated(gjorb3%i1)) call mma_deallocate(gjorb3%i1)
  if (allocated(gjorb3%i2)) call mma_deallocate(gjorb3%i2)
  if (allocated(cvbstot)) call mma_deallocate(cvbstot)
  if (allocated(cvbsspn)) call mma_deallocate(cvbsspn)
  if (allocated(cvbdet)) call mma_deallocate(cvbdet)
  if (allocated(dvbdet)) call mma_deallocate(dvbdet)
  if (allocated(evbdet)) call mma_deallocate(evbdet)
  call mma_deallocate(orbstry)
  call mma_deallocate(cvbtry)
  call mma_deallocate(vbdet)
end if
release(7) = .true.
release(8) = .false.

if ((icase7 == 1) .or. (icase7 == 2) .or. (icase7 == 3)) then
  !FIXME: These deallocations should not be necessary
  if (allocated(orbinv)) call mma_deallocate(orbinv)
  if (allocated(sorbs)) call mma_deallocate(sorbs)
  if (allocated(owrk2)) call mma_deallocate(owrk2)
  if (allocated(gjorb%r)) call mma_deallocate(gjorb%r)
  if (allocated(gjorb%i1)) call mma_deallocate(gjorb%i1)
  if (allocated(gjorb%i2)) call mma_deallocate(gjorb%i2)
  if (allocated(gjorb2%r)) call mma_deallocate(gjorb2%r)
  if (allocated(gjorb2%i1)) call mma_deallocate(gjorb2%i1)
  if (allocated(gjorb2%i2)) call mma_deallocate(gjorb2%i2)
  if (allocated(gjorb3%r)) call mma_deallocate(gjorb3%r)
  if (allocated(gjorb3%i1)) call mma_deallocate(gjorb3%i1)
  if (allocated(gjorb3%i2)) call mma_deallocate(gjorb3%i2)
  if (allocated(cvbstot)) call mma_deallocate(cvbstot)
  if (allocated(cvbsspn)) call mma_deallocate(cvbsspn)
  if (allocated(cvbdet)) call mma_deallocate(cvbdet)
  if (allocated(dvbdet)) call mma_deallocate(dvbdet)
  if (allocated(evbdet)) call mma_deallocate(evbdet)
  call mma_allocate(orbinv,norb,norb,label='orbinv')
  call mma_allocate(sorbs,norb,norb,label='sorbs')
  call mma_allocate(owrk2,norb,norb,label='owrk2')
  call mma_allocate(gjorb%r,norb,norb,label='gjorb%r')
  call mma_allocate(gjorb%i1,norb,label='gjorb%i1')
  call mma_allocate(gjorb%i2,2,norb*norb,label='gjorb%i2')
  call mma_allocate(gjorb2%r,norb,norb,label='gjorb2%r')
  call mma_allocate(gjorb2%i1,norb,label='gjorb2%i1')
  call mma_allocate(gjorb2%i2,2,norb*norb,label='gjorb2%i2')
  call mma_allocate(gjorb3%r,norb,norb,label='gjorb3%r')
  call mma_allocate(gjorb3%i1,norb,label='gjorb3%i1')
  call mma_allocate(gjorb3%i2,2,norb*norb,label='gjorb3%i2')
  call mma_allocate(cvbstot,nvb,label='cvbstot')
  call mma_allocate(cvbsspn,nvb,label='cvbsspn')
  call mma_allocate(cvbdet,ndetvb,label='cvbdet')
  call mma_allocate(dvbdet,ndetvb,label='dvbdet')
  call mma_allocate(evbdet,ndetvb,label='evbdet')
end if

!FIXME: These deallocations should not be necessary
if (allocated(orbstry)) call mma_deallocate(orbstry)
if (allocated(cvbtry)) call mma_deallocate(cvbtry)
if (allocated(vbdet)) call mma_deallocate(vbdet)
call mma_allocate(orbstry,norb,norb,label='orbstry')
call mma_allocate(cvbtry,ndetvb,label='cvbtry')
call mma_allocate(vbdet,ndetvb,label='vbdet')

return

end subroutine chop7_cvb
