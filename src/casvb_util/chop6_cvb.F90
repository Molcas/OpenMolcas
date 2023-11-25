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

use casvb_global, only: endvar, grad1, grad2, gradx, hessorb, hesst, icase6, icnt_ci, ishstruc, ivbweights, maxdav, mxdav, &
                        mxirrep, mxiter, ncivb, ndet, ndetvb, nirrep, norb, npr, nprorb, nvb, release, sstruc, sstruc2, variat, &
                        vec1, wdx
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6, RtoI

implicit none
integer(kind=iwp) :: idav, ir, iremain, mem_applyh, memwrk, ncimx, need
logical(kind=iwp) :: done

if (release(6)) then
  if (allocated(sstruc)) call mma_deallocate(sstruc)
  if (allocated(sstruc2)) call mma_deallocate(sstruc2)
  if (allocated(hessorb)) call mma_deallocate(hessorb)
  if (allocated(hesst)) call mma_deallocate(hesst)
  if (allocated(wdx)) call mma_deallocate(wdx)
  if (allocated(grad1)) call mma_deallocate(grad1)
  if (allocated(grad2)) call mma_deallocate(grad2)
  if (allocated(gradx)) call mma_deallocate(gradx)
  if (allocated(vec1)) call mma_deallocate(vec1)
end if
release(6) = .true.
release(7) = .false.

icnt_ci(6) = 0
if (icase6 == 1) then
  !FIXME: These deallocations should not be needed
  if (allocated(sstruc)) call mma_deallocate(sstruc)
  if (allocated(sstruc2)) call mma_deallocate(sstruc2)
  if (allocated(hessorb)) call mma_deallocate(hessorb)
  if (allocated(hesst)) call mma_deallocate(hesst)
  if (allocated(wdx)) call mma_deallocate(wdx)
  if (allocated(grad1)) call mma_deallocate(grad1)
  if (allocated(grad2)) call mma_deallocate(grad2)
  if (allocated(gradx)) call mma_deallocate(gradx)
  if (allocated(vec1)) call mma_deallocate(vec1)
  ! Standard non-linear optimization procedure:
  call mma_allocate(sstruc,norb*norb+nvb+1+mxirrep,1,label='sstruc')
  call mma_allocate(sstruc2,npr,1,label='sstruc2')
  call mma_allocate(hessorb,nprorb,nprorb,label='hessorb')
  call mma_allocate(hesst,norb**2,norb**2,label='hesst')
  call mma_allocate(wdx,npr,label='wdx')
  call mma_allocate(grad1,npr,label='grad1')
  call mma_allocate(grad2,npr,label='grad2')
  call mma_allocate(gradx,norb,norb,label='gradx')
  call mma_allocate(vec1,max(npr,ndetvb),label='vec1')
else if (icase6 == 2) then
  ! Overlap-based Davidson optimization:
  call mma_maxDBLE(iremain)
  maxdav = min(mxiter,nvb,mxdav)

  memwrk = ndetvb+5*norb*norb+3*(norb+2*norb*norb+RtoI-1)/RtoI
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
    write(u6,*) ' Not enough memory for Davidson!',need,iremain
    call abend_cvb()
  end if
  maxdav = idav

else if (icase6 == 3) then
  ! Energy-based Davidson optimization:
  call mma_maxDBLE(iremain)
  maxdav = min(mxiter,nvb,mxdav)

  ! neread is uninitialized
  !mem_applyh = ndet+neread
  mem_applyh = ndet
  ncimx = 0
  do ir=1,nirrep
    ncimx = max(ncimx,ncivb(ir))
  end do
  if (ncimx /= ndet) mem_applyh = mem_applyh+ncimx
  memwrk = ndetvb+3*norb*norb+2*(norb+2*norb*norb+RtoI-1)/RtoI

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
    write(u6,*) ' Not enough memory for Davidson!',need,iremain
    call abend_cvb()
  end if
  maxdav = idav

else if (icase6 == 4) then
  ! Wavefunction analysis:
  if (((.not. variat) .or. endvar) .and. ((ivbweights > 1) .or. (ishstruc == 1))) then
    !FIXME: These deallocations should not be needed
    if (allocated(sstruc)) call mma_deallocate(sstruc)
    if (allocated(sstruc2)) call mma_deallocate(sstruc2)
    call mma_allocate(sstruc,nvb,nvb,label='sstruc')
    call mma_allocate(sstruc2,nvb,nvb,label='sstruc2')
  end if
end if

return

end subroutine chop6_cvb
