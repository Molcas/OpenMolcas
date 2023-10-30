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

subroutine symgen_cvb(nalf1,nbet1,nda1,ndb1,isymalf,isymbet,iasyind,ibsyind,irpdet)

use Symmetry_Info, only: Mul
use casvb_global, only: ityp, mxirrep, norb
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nalf1, nbet1, nda1, ndb1
integer(kind=iwp), intent(out) :: isymalf(nda1), isymbet(ndb1), iasyind(0:mxirrep), ibsyind(0:mxirrep), irpdet(mxirrep)
integer(kind=iwp) :: ia, ib, icount(mxirrep), ida, idb, indx, iorb, irp, irpalf(mxirrep), irpbet(mxirrep), irrep, jrp, rc
integer(kind=iwp), allocatable :: ialfsym(:), ibetsym(:), locc(:), lunocc(:), maxgrph(:), mingrph(:), nk(:), xalf(:,:), xbet(:,:)

call mma_allocate(mingrph,[0,norb],label='mingrph')
call mma_allocate(maxgrph,[0,norb],label='maxgrph')
call mma_allocate(nk,[0,norb],label='nk')
call mma_allocate(locc,norb+1,label='locc')
call mma_allocate(lunocc,norb+1,label='lunocc')

! Alpha loop:
call mma_allocate(ialfsym,nda1,label='ialfsym')
call mma_allocate(xalf,[0,norb],[0,nalf1],label='xalf')
irpalf(:) = 0
do iorb=0,norb
  mingrph(iorb) = max(iorb-norb+nalf1,0)
  maxgrph(iorb) = min(iorb,nalf1)
end do
call weight_cvb(xalf,mingrph,maxgrph,nalf1,norb)
nk(:) = maxgrph(:)
call occupy_cvb(nk,norb,locc,lunocc)
indx = 1
do
  irp = 1
  do ia=1,nalf1
    irp = Mul(irp,ityp(locc(ia)))
  end do
  irpalf(irp) = irpalf(irp)+1
  ialfsym(indx) = irp
  call loind_cvb(norb,nalf1,nk,mingrph,maxgrph,locc,lunocc,indx,xalf,rc)
  if (rc == 0) exit
end do
iasyind(0) = 0
do irp=1,mxirrep
  iasyind(irp) = iasyind(irp-1)+irpalf(irp)
end do
icount(:) = 0
do ida=1,nda1
  irrep = ialfsym(ida)
  icount(irrep) = icount(irrep)+1
  isymalf(icount(irrep)+iasyind(irrep-1)) = ida
end do
call mma_deallocate(ialfsym)
call mma_deallocate(xalf)

! Beta loop:
call mma_allocate(ibetsym,ndb1,label='ibetsym')
call mma_allocate(xbet,[0,norb],[0,nbet1],label='xbet')
irpbet(:) = 0
do iorb=0,norb
  mingrph(iorb) = max(iorb-norb+nbet1,0)
  maxgrph(iorb) = min(iorb,nbet1)
end do
call weight_cvb(xbet,mingrph,maxgrph,nbet1,norb)
nk(:) = maxgrph(:)
call occupy_cvb(nk,norb,locc,lunocc)
indx = 1
do
  irp = 1
  do ib=1,nbet1
    irp = Mul(irp,ityp(locc(ib)))
  end do
  irpbet(irp) = irpbet(irp)+1
  ibetsym(indx) = irp
  call loind_cvb(norb,nbet1,nk,mingrph,maxgrph,locc,lunocc,indx,xbet,rc)
  if (rc == 0) exit
end do
ibsyind(0) = 0
do irp=1,mxirrep
  ibsyind(irp) = ibsyind(irp-1)+irpbet(irp)
end do
icount(:) = 0
do idb=1,ndb1
  irrep = ibetsym(idb)
  icount(irrep) = icount(irrep)+1
  isymbet(icount(irrep)+ibsyind(irrep-1)) = idb
end do
call mma_deallocate(ibetsym)
call mma_deallocate(xbet)

call mma_deallocate(mingrph)
call mma_deallocate(maxgrph)
call mma_deallocate(nk)
call mma_deallocate(locc)
call mma_deallocate(lunocc)

irpdet(:) = 0
do irp=1,mxirrep
  do jrp=1,mxirrep
    irpdet(irp) = irpdet(irp)+irpalf(jrp)*irpbet(Mul(irp,jrp))
  end do
end do

return

end subroutine symgen_cvb
