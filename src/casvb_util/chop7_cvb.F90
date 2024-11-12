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
  call mma_deallocate(orbinv,safe='*')
  call mma_deallocate(sorbs,safe='*')
  call mma_deallocate(owrk2,safe='*')
  call mma_deallocate(gjorb%r,safe='*')
  call mma_deallocate(gjorb%i1,safe='*')
  call mma_deallocate(gjorb%i2,safe='*')
  call mma_deallocate(gjorb2%r,safe='*')
  call mma_deallocate(gjorb2%i1,safe='*')
  call mma_deallocate(gjorb2%i2,safe='*')
  call mma_deallocate(gjorb3%r,safe='*')
  call mma_deallocate(gjorb3%i1,safe='*')
  call mma_deallocate(gjorb3%i2,safe='*')
  call mma_deallocate(cvbstot,safe='*')
  call mma_deallocate(cvbsspn,safe='*')
  call mma_deallocate(cvbdet,safe='*')
  call mma_deallocate(dvbdet,safe='*')
  call mma_deallocate(evbdet,safe='*')
  call mma_deallocate(orbstry)
  call mma_deallocate(cvbtry)
  call mma_deallocate(vbdet)
end if
release(7) = .true.
release(8) = .false.

if ((icase7 == 1) .or. (icase7 == 2) .or. (icase7 == 3)) then
  !FIXME: These deallocations should not be necessary
  call mma_deallocate(orbinv,safe='*')
  call mma_deallocate(sorbs,safe='*')
  call mma_deallocate(owrk2,safe='*')
  call mma_deallocate(gjorb%r,safe='*')
  call mma_deallocate(gjorb%i1,safe='*')
  call mma_deallocate(gjorb%i2,safe='*')
  call mma_deallocate(gjorb2%r,safe='*')
  call mma_deallocate(gjorb2%i1,safe='*')
  call mma_deallocate(gjorb2%i2,safe='*')
  call mma_deallocate(gjorb3%r,safe='*')
  call mma_deallocate(gjorb3%i1,safe='*')
  call mma_deallocate(gjorb3%i2,safe='*')
  call mma_deallocate(cvbstot,safe='*')
  call mma_deallocate(cvbsspn,safe='*')
  call mma_deallocate(cvbdet,safe='*')
  call mma_deallocate(dvbdet,safe='*')
  call mma_deallocate(evbdet,safe='*')
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
call mma_deallocate(orbstry,safe='*')
call mma_deallocate(cvbtry,safe='*')
call mma_deallocate(vbdet,safe='*')
call mma_allocate(orbstry,norb,norb,label='orbstry')
call mma_allocate(cvbtry,ndetvb,label='cvbtry')
call mma_allocate(vbdet,ndetvb,label='vbdet')

return

end subroutine chop7_cvb
