!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine cipro()

use gugaci_global, only: denm1, LuCiDen, LuCiMO, max_root, maxgdm, mroot, ng_sm, nlsm_all, nlsm_bas, pror
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), parameter :: maxmolcasorb = 5000, maxpro = 50
integer(kind=iwp) :: i, icall, idisk, idummy(1), idx_idisk0(64), idx_idisk1(max_root+1), iend, im, iopt, ipc, ipcom(maxpro), &
                     iprop, irec, iroot, irtc, ista, isymlb, nc, nc0, nc1, nc2, nlsm_del(maxgdm), nmo, npro, nsiz
real(kind=wp) :: pgauge(3,maxpro), pnuc(maxpro)
character :: bsbl(2*4*maxmolcasorb)
character(len=8) :: label, pname(maxpro), ptyp(maxpro)
real(kind=wp), allocatable :: cmo(:), cno(:), occ(:)
real(kind=wp), pointer :: denao(:), omat(:), vprop(:,:,:)

! mrci nature orbital
idisk = 0
call idafile(lucimo,2,idx_idisk0,64,idisk)
!write(u6,*) 'read head'
idisk = idx_idisk0(2)
nc = 2*4*maxmolcasorb
call cdafile(lucimo,2,bsbl,nc,idisk)

!open(100,file='dat1')
nc0 = 0
nc1 = 0
nc2 = 0
nmo = 0
do im=1,ng_sm
  nc0 = nc0+nlsm_bas(im)**2
  nc1 = nc1+nlsm_all(im)*(nlsm_all(im)+1)/2
  nc2 = nc2+nlsm_bas(im)*(nlsm_bas(im)+1)/2
  nmo = nmo+nlsm_bas(im)
end do
! read property labels
iopt = 8
npro = 0
do i=1,100
  label = 'undef'
  call irdone(irec,iopt+1,label,ipc,idummy,isymlb)
  if (irec /= 0) exit
  iopt = 16
  if (mod(isymlb,2) == 0) cycle
  npro = npro+1
  pname(npro) = label
  ipcom(npro) = ipc
  ptyp(npro) = 'HERM'
  if (label == 'VELOCITY') ptyp(npro) = 'ANTI'
  if (label == 'ANGMOM  ') ptyp(npro) = 'ANTI'
  if (label(1:5) == 'MLTPV') ptyp(npro) = 'ANTI'
  if (npro >= maxpro) exit
end do

!write(u6,'(10(1X,a8))') pname(1:npro)
!write(u6,*) ipcom(1:npro)
!write(u6,*) ptyp(1:npro)

allocate(omat(nc2))
allocate(denao(nc0))
allocate(vprop(mroot,mroot,npro))
!allocate(denao(nmo,nmo))
call mma_allocate(cmo,nc0,label='cmo')
call mma_allocate(cno,nc0,label='cno')
call mma_allocate(occ,nmo,label='occ')
idisk = idx_idisk0(3)
call ddafile(lucimo,2,cmo,nc0,idisk)
! read overlap matrix
call rdone(irec,6,'MLTPL  0',1,omat,idummy(1))
idisk = 0
call idafile(luciden,2,idx_idisk1,max_root+1,idisk)

icall = 0
nlsm_del(:) = 0
do iroot=1,mroot
  ! read density matrix
  idisk = idx_idisk1(iroot)
  call ddafile(luciden,2,denm1,nc1,idisk)
  call natureorb(nlsm_bas,nlsm_all,nlsm_del,ng_sm,denm1,nc1,cmo,nc0,bsbl,nc,cno,occ,nmo,pror)

  ! mulliken population
  call xflush(u6)
  write(u6,'(a,i2)') ' mulliken charges for state nr ',iroot
  !call charge(nsym,nbas,name,cno,occ,smat,2,.true.,.true.)
  call charge(ng_sm,nlsm_bas,bsbl,cno,occ,omat,2,.true.,.false.)
  write(u6,*) ' ',('*',i=1,70)
  call xflush(u6)
  ! transform mo density matrix to ao basis
  !write(u6,'(10i8)') nc0,nmo,nlsm_bas
  call transden(ng_sm,nlsm_bas,denao,cno,nc0,occ,nmo)

  ! calculated properties
  call calprop(ng_sm,nlsm_bas,mroot,iroot,iroot,nc2,npro,pname,ipcom,ptyp,denao,nc0,vprop,pgauge,pnuc,icall)
  ! print property
end do
call mma_deallocate(cmo)
call mma_deallocate(cno)
call mma_deallocate(occ)
!close(100)

! this code copy from molcas
if (npro > 0) then
  ! write expectation values:
  write(u6,*)
  write(u6,*) ' expectation values of various operators:'
  write(u6,*) '(note: electronic multipoles include a negative sign.)'
  do iprop=1,npro
    if (ptyp(iprop) == 'ANTI') cycle
    do ista=1,mroot,4
      iend = min(ista+3,mroot)
      write(u6,*)
      write(u6,'(1x,a,a8,a,i4)') '   property: ',pname(iprop),'   component:',ipcom(iprop)
      write(u6,'(1x,a,3f16.8)') '    gauge origin:',(pgauge(i,iprop),i=1,3)
      write(u6,'(1x,a,i8,3i16)') '            root:',(i,i=ista,iend)
      write(u6,'(1x,a,4f16.8)') '      electronic:',(vprop(i,i,iprop),i=ista,iend)
      write(u6,'(1x,a,4f16.8)') '         nuclear:',(pnuc(iprop),i=ista,iend)
      write(u6,'(1x,a,4f16.8)') '           total:',(pnuc(iprop)+vprop(i,i,iprop),i=ista,iend)
    end do
  end do
  write(u6,*)
end if

iopt = 8
npro = 0
do i=1,100
  label = 'undef'
  call irdone(irec,iopt+1,label,ipc,idummy,isymlb)
  if (irec /= 0) exit
  iopt = 16
  !if (mod(isymlb,2) == 0) cycle
  npro = npro+1
  pname(npro) = label
  ipcom(npro) = ipc
  ptyp(npro) = 'HERM'
  if (label == 'VELOCITY') ptyp(npro) = 'ANTI'
  if (label == 'ANGMOM  ') ptyp(npro) = 'ANTI'
  if (label(1:5) == 'MLTPV') ptyp(npro) = 'ANTI'
  if (npro >= maxpro) exit
end do

!write(u6,'(10(1x,a8))') pname(1:npro)
!write(u6,'(10i4)') ipcom(1:npro)
!write(u6,'(10(1x,a8))') ptyp(1:npro)

nsiz = 0
call Molcas_BinaryOpen_Vanilla(110,'soint.dat')
do i=1,npro
  if (pname(i)(1:4) /= 'AMFI') cycle
  omat = Zero
  call irdone(irtc,1,pname(i),ipcom(i),idummy,isymlb)
  if (irtc == 0) nsiz = idummy(1)
  call rdone(irtc,0,pname(i),ipcom(i),omat,isymlb)
  if (nsiz > nc2) then
    write(u6,*) 'in subroutine cipro, read so int error'
    call abend()
  end if
  write(u6,'(a8,1x,2i4)') 'nsiz=',i,nsiz
  write(u6,'(5(1x,f12.8))') omat(1:nsiz)
  write(110) nsiz
  write(110) omat(1:nsiz)
end do
close(110)
deallocate(omat)
deallocate(denao)
deallocate(vprop)

return

end subroutine cipro

subroutine calprop(ngsm,nsbas,mroot,istate,jstate,nsi,npro,pname,ipcom,ptyp,aden,lmo,vprop,pgauge,pnuc,icall)

use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6, r8

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ngsm, nsbas(ngsm), mroot, istate, jstate, nsi, npro, ipcom(npro), lmo
character(len=8), intent(_IN_) :: pname(npro), ptyp(npro)
real(kind=wp), intent(in) :: aden(lmo)
real(kind=wp), intent(inout) :: vprop(mroot,mroot,npro)
real(kind=wp), intent(inout) :: pgauge(3,npro), pnuc(npro)
integer(kind=iwp), intent(inout) :: icall
integer(kind=iwp) :: i, idummy(1), im, irtc, isymlb, j, nc, nc0, nc1, nsiz
real(kind=wp) :: amat(nsi), pint(nsi+4), sgn, smat(nsi), val
real(kind=r8), external :: ddot_

! we have two kind of density matrix, symm or anti-symm
! compress density matrix
nc0 = 0
nc1 = 1
nc = 0
do im=1,ngsm
  if (nsbas(im) == 0) cycle
  do i=1,nsbas(im)
    do j=1,i-1
      smat(nc1) = aden(nc0+(j-1)*nsbas(im)+i)+aden(nc0+(i-1)*nsbas(im)+j)
      amat(nc1) = aden(nc0+(j-1)*nsbas(im)+i)-aden(nc0+(i-1)*nsbas(im)+j)
      !smat(nc1) = aden(nc+i,nc+j)+aden(nc+j,nc+i)
      !amat(nc1) = aden(nc+i,nc+j)-aden(nc+j,nc+i)
      nc1 = nc1+1
    end do
    smat(nc1) = aden(nc0+(i-1)*nsbas(im)+i)
    !smat(nc1) = aden(nc+i,nc+i)
    amat(nc1) = Zero
    nc1 = nc1+1
  end do
  nc = nc+nsbas(im)
  nc0 = nc0+nsbas(im)**2
end do

!open(100,file='dat1')
!do i=1,1293
!  write(100,'(i8,1x,f14.8)') i,smat(i)
!end do
!close(100)

nsiz = nsi
! read property int and calculated property
do i=1,npro
  call irdone(irtc,1,pname(i),ipcom(i),idummy,isymlb)
  if (irtc == 0) nsiz = idummy(1)
  call rdone(irtc,0,pname(i),ipcom(i),pint,isymlb)
  !write(u6,*) 'nsiz',nsiz
  if (icall == 0) then
    pgauge(1,i) = pint(nsiz+1)
    pgauge(2,i) = pint(nsiz+2)
    pgauge(3,i) = pint(nsiz+3)
    pnuc(i) = pint(nsiz+4)
  end if
  if (isymlb /= 1) then
    write(u6,*) 'error calcualte property,need debug'
    call abend()
  end if
  !write(u6,*) 'pop int'
  !write(u6,'(10(1x,f8.5))') pint(1:100)
  !write(u6,*) 'sden '
  !write(u6,'(10(1x,f8.5))') smat(1:100)

  ! calculate and save property
  sgn = One
  if (pname(i)(1:5) == 'MLTPL') sgn = -One
  if (ptyp(i)(1:4) == 'HERM') then
    val = sgn*ddot_(nsiz,smat,1,pint,1)
    vprop(istate,jstate,i) = val
    vprop(jstate,istate,i) = val
  else
    val = sgn*ddot_(nsiz,amat,1,pint,1)
    vprop(istate,jstate,i) = val
    vprop(jstate,istate,i) = -val
  end if
  !write(u6,*) nsiz,'val=',val
end do

icall = 1

return

end subroutine calprop

subroutine transden(ngsm,nsbas,denao,cno,lmo,occ,loc)
! calculate density matrix in ao basis

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ngsm, nsbas(ngsm), lmo, loc
real(kind=wp), intent(out) :: denao(lmo)
real(kind=wp), intent(in) :: cno(lmo), occ(loc)
integer(kind=iwp) :: i, im, nc, nc0, nc1, ni
real(kind=wp) :: val

denao(:) = Zero
!write(u6,*) 'occ',nsbas(1:ngsm)
!write(u6,'(10(1x,f8.5))') occ
nc = 1
nc0 = 0
nc1 = 1
do im=1,ngsm
  ni = nsbas(im)
  if (ni == 0) cycle
  !write(u6,*) nc,ni
  do i=1,ni
    val = occ(i+nc0)
    call dger(ni,ni,val,cno(nc1),1,cno(nc1),1,denao(nc),ni)
    nc1 = nc1+ni
  end do
  nc0 = nc0+ni
  nc = nc+ni**2
end do

return

end subroutine transden
