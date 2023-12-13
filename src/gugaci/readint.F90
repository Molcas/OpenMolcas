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
! Copyright (C) 2007, Bingbing Suo                                     *
!***********************************************************************
!  16 apr 2007 - bsuo - revised to use molcas integrals

#ifdef MOLPRO

subroutine intrd()

use gugaci_global, only: lsmorb, LuTwoMO, map_orb_order, nlsm_all, nlsm_bas, ng_sm, noidx, voint, vpotnuc
use file_qininit, only: maxrecord
use Symmetry_Info, only: Mul
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: i, idx, lrcii, lrcij, lri, lrj, lrt, nc, ni, nidx, nintone, nism, nmob, nsmint
real(kind=wp) :: ecor
integer(kind=iwp), allocatable :: noffset(:)
real(kind=wp), allocatable :: xfock(:)

nintone = 0
nmob = 0
nidx = 0
ni = 0
do i=1,ng_sm
  nism = nlsm_all(i)
  nsmint = nism*(nism+1)/2
  nmob = nmob+nlsm_bas(i)*nlsm_bas(i)
  noidx(i) = nidx
  nidx = nidx+nism
  nintone = nintone+nsmint
  lsmorb(ni+1:nidx) = i
  ni = nidx
end do

call mma_allocate(noffset,maxrecord,label='noffset')
noffset(:) = 0
ni = 0
call idafile(Lutwomo,2,noffset,maxrecord,ni)
! 2nd record in file traint, nuclear repulsive energy
call ddafile(Lutwomo,2,vpotnuc,1,ni)
! 3rd record, ecore
call ddafile(Lutwomo,2,ecor,1,ni)  ! effective core energy
write(u6,'(a,1x,f14.8)') ' Nuclear repulsive energy:',vpotnuc
write(u6,'(a,1x,f14.8)') ' Frozen  core energy:     ',ecor
vpotnuc = vpotnuc+ecor
! 4th record, 1e int
call mma_allocate(xfock,nintone,label='xfock')
call ddafile(Lutwomo,2,xfock,nintone,ni) ! 1e integrals

! write one electron fock matrix into voint
nidx = 0
do i=1,ng_sm
  nism = nlsm_all(i)
  idx = noidx(i)
  nsmint = nism*(nism+1)/2
  nc = 0
  do lri=1,nism
    do lrj=1,lri
      nc = nc+1
      lrcii = map_orb_order(lri+idx)
      lrcij = map_orb_order(lrj+idx)
      if (lrcii > lrcij) then
        lrt = lrcij
        lrcij = lrcii
        lrcii = lrt
      end if
      voint(lrcii,lrcij) = xfock(nc+nidx)
      !write(u6,'(1x,i3,1x,i3,1x,f18.9)') lrcii,lrcij,xfock(nc+nidx)
    end do
  end do
  nidx = nidx+nsmint
end do
call readtwoeint_molpro(lutwomo,maxrecord,noffset,nlsm_all,ng_sm,Mul,map_orb_order,noidx)
call mma_deallocate(noffset)
call mma_deallocate(xfock)

!write(u6,*)
!write(u6,*) 'MRCI integrals'
!do lri=1,406
!  write(55,'(i8,f16.8)') lri,vector1(lri)
!end do

return

end subroutine intrd

subroutine readtwoeint_molpro(nft,maxrecord,noffset,norb,ngsm,multab,maporb,noidx)

use gugaci_global, only: max_orb, ntrabuf
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nft, maxrecord, noffset(maxrecord), norb(8), ngsm, multab(8,8), maporb(max_orb), noidx(8)
integer(kind=iwp) :: idisk, idx, iout, ityp, li, lj, lk, ll, lri, lrj, lrk, lrl, nbpq, nbrs, nintb, nop, noq, nor, nos, nsp, nspq, &
                     nspqr, nsq, nsr, nss, nssm, ntj, ntk, ntl
real(kind=wp) :: val
real(kind=wp), allocatable :: buff(:)

idisk = noffset(5)
write(u6,2000)
call mma_allocate(buff,ntrabuf,label='buff')
do nsp=1,ngsm
  nop = norb(nsp)
  do nsq=1,nsp
    noq = norb(nsq)
    nspq = multab(nsp,nsq)
    do nsr=1,nsp
      nor = norb(nsr)
      nspqr = multab(nspq,nsr)
      nssm = nsr
      if (nsr == nsp) nssm = nsq
      do nss=1,nssm
        if (nspqr /= nss) cycle
        nos = norb(nss)

        ityp = 0
        if (nsr == nss) then
          nbpq = (nop+nop**2)/2
          nbrs = (nos+nos**2)/2
          if (nsp == nsr) then
            ! (ii|ii) type 1 int
            nintb = (nbpq+nbpq**2)/2
            ityp = 1
          else
            ! (ii|jj) type 3 int
            nintb = nbpq*nbrs
            ityp = 3
          end if
        else
          nbpq = nop*noq
          nbrs = nor*nos
          if (nsp == nsr) then
            ! (ij|ij) type 2 int
            nintb = (nbpq+nbpq**2)/2
            ityp = 2
          else
            ! (ij|kl) type 4 int
            nintb = nbpq*nbrs
            ityp = 4
          end if
        end if

        if (nintb == 0) cycle
        write(u6,2100) nsp,nsq,nsr,nss,nop,noq,nor,nos,nintb

        idx = 0
        iout = 0
        !call ddatard(nft,buff,ntrabuf,idisk)
        if (nintb > ntrabuf) then
          call ddafile(nft,2,buff,ntrabuf,idisk)
        else
          call ddafile(nft,2,buff,nintb,idisk)
        end if

        do li=1,nos
          ntj = nor
          if ((ityp == 1) .or. (ityp == 3)) ntj = li
          do lj=1,ntj
            ntk = noq
            if ((ityp == 1) .or. (ityp == 2)) ntk = li
            do lk=1,ntk
              ntl = nop
              if ((ityp == 1) .or. (ityp == 3)) ntl = lk
              if ((ityp == 1) .and. (li == lk)) ntl = lj
              if ((ityp == 2) .and. (li == lk)) ntl = lj
              do ll=1,ntl
                iout = iout+1
                if (iout > ntrabuf) then
                  if (nintb-idx < ntrabuf) then
                    call ddafile(nft,2,buff,nintb-idx,idisk)
                  else
                    call ddafile(nft,2,buff,ntrabuf,idisk)
                  end if
                  iout = 1
                end if
                idx = idx+1
                val = buff(iout)
                lri = maporb(li+noidx(nss))
                lrj = maporb(lj+noidx(nsr))
                lrk = maporb(lk+noidx(nsq))
                lrl = maporb(ll+noidx(nsp))
                !if (ityp /= 2) cycle
                !write(u6,'(9(1x,i2),1x,f18.9)') li,lj,lk,ll,lri,lrj,lrk,lrl,iout,val
                call intrw_mol(lri,lrj,lrk,lrl,val)
              end do
            end do
          end do
        end do

      end do
    end do
  end do
end do
call mma_deallocate(buff)

return

2000 format(/7x,'symmetry',6x,' orbitals',8x,'integrals')
2100 format(7x,4i2,1x,4i4,2x,3x,i9)

end subroutine readtwoeint_molpro

#else

subroutine intrd_molcas()

use gugaci_global, only: FnOneMO, FnTwoMO, lsmorb, LuOneMO, LuTwoMO, map_orb_order, ng_sm, nlsm_all, nlsm_bas, noidx, voint, vpotnuc
use Symmetry_Info, only: Mul
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: i, idisk, idx, lrcii, lrcij, lri, lrj, lrt, nc, ni, nidx, nintone, nism, nmob, nsmint
real(kind=wp) :: ecor
real(kind=wp), allocatable :: x(:), xfock(:)

nintone = 0
nmob = 0
nidx = 0
ni = 0
do i=1,ng_sm
  nism = nlsm_all(i)
  nsmint = nism*(nism+1)/2
  nmob = nmob+nlsm_bas(i)*nlsm_bas(i)
  noidx(i) = nidx
  nidx = nidx+nism
  nintone = nintone+nsmint
  lsmorb(ni+1:nidx) = i
  ni = nidx
end do
call mma_allocate(x,nmob,label='x')
call daname(luonemo,fnonemo)
call readtraonehead(luonemo,ecor,idisk)
vpotnuc = ecor

! read mo coeff, need debug, if frozen and delete orbitals are not zero
!call ddatard(nft,x,nmob,idisk)
call ddafile(luonemo,2,x,nmob,idisk)

! read one electron fock matrix
call mma_allocate(xfock,nintone,label='xfock')
call ddafile(luonemo,2,xfock,nintone,idisk)
!call ddatard(nft,xfock,nintone,idisk)
call daclos(luonemo)
! read one electron kinetic integrals
!call ddatard(nft,x1e,nintone,idisk)
!call cclose_molcas(nft)
! write one electron fock matrix into voint
nidx = 0
do i=1,ng_sm
  nism = nlsm_all(i)
  idx = noidx(i)
  nsmint = nism*(nism+1)/2
  nc = 0
  do lri=1,nism
    do lrj=1,lri
      nc = nc+1
      lrcii = map_orb_order(lri+idx)
      lrcij = map_orb_order(lrj+idx)
      if (lrcii > lrcij) then
        lrt = lrcij
        lrcij = lrcii
        lrcii = lrt
      end if
      voint(lrcii,lrcij) = xfock(nc+nidx)
      !write(u6,'(1x,i3,1x,i3,1x,f18.9)') lrcii,lrcij,xfock(nc+nidx)
    end do
  end do
  nidx = nidx+nsmint
end do
call mma_deallocate(xfock)

call daname(lutwomo,fntwomo)
call readtwoeint(lutwomo,nlsm_all,ng_sm,Mul,map_orb_order,noidx)
call daclos(lutwomo)
write(u6,*)
call mma_deallocate(x)

return

end subroutine intrd_molcas

subroutine readtwoeint(nft,norb,ngsm,multab,maporb,noidx)

use gugaci_global, only: lenintegral, max_orb, ntrabuf, ntratoc
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nft, norb(8), ngsm, multab(8,8), maporb(max_orb), noidx(8)
integer(kind=iwp) :: idisk, idx, iout, lenrd, li, lj, lk, ll, lri, lrj, lrk, lrl, nbpq, nbrs, nintb, nop, noq, nor, nos, nsp, &
                     nspq, nspqr, nsq, nsr, nss, nssm, ntj, ntk, numax, numin
real(kind=wp) :: val
integer(kind=iwp), allocatable :: itratoc(:)
real(kind=wp), allocatable :: buff(:)

idisk = 0
lenrd = ntratoc*lenintegral
write(u6,*) lenrd
call mma_allocate(itratoc,ntratoc,label='itratoc')
!call idatard(nft,itratoc,lenrd,idisk)
call idafile(nft,2,itratoc,ntratoc,idisk)
!write(u6,'(10(1x,i8))') itratoc(1:10)
call mma_deallocate(itratoc)
call mma_allocate(buff,ntrabuf,label='buff')
write(u6,2000)
do nsp=1,ngsm
  nop = norb(nsp)
  do nsq=1,nsp
    noq = norb(nsq)
    nspq = multab(nsp,nsq)
    do nsr=1,nsp
      nor = norb(nsr)
      nspqr = multab(nspq,nsr)
      nssm = nsr
      if (nsr == nsp) nssm = nsq
      do nss=1,nssm
        if (nspqr /= nss) cycle
        nos = norb(nss)

        if (nsr == nss) then
          nbpq = (nop+nop**2)/2
          nbrs = (nos+nos**2)/2
          if (nsp == nsr) then
            ! (ii|ii) type 1 int
            nintb = (nbpq+nbpq**2)/2
          else
            ! (ii|jj) type 3 int
            nintb = nbpq*nbrs
          end if
        else
          nbpq = nop*noq
          nbrs = nor*nos
          if (nsp == nsr) then
            ! (ij|ij) type 2 int
            nintb = (nbpq+nbpq**2)/2
          else
            ! (ij|kl) type 4 int
            nintb = nbpq*nbrs
          end if
        end if

        if (nintb == 0) cycle
        write(u6,2100) nsp,nsq,nsr,nss,nop,noq,nor,nos,nintb

        idx = 0
        iout = 0
        !call ddatard(nft,buff,ntrabuf,idisk)
        call ddafile(nft,2,buff,ntrabuf,idisk)

        do li=1,nor
          ntj = nos
          if (nsr == nss) ntj = li
          do lj=1,ntj
            ntk = 1
            if (nsp == nsr) ntk = li
            do lk=ntk,nop
              numin = 1
              if ((nsp == nsr) .and. (lk == li)) numin = lj
              numax = noq
              if (nsp == nsq) numax = lk
              do ll=numin,numax
                iout = iout+1
                if (iout > ntrabuf) then
                  !call ddatard(nft,buff,ntrabuf,idisk)
                  call ddafile(nft,2,buff,ntrabuf,idisk)
                  iout = 1
                end if
                idx = idx+1
                val = buff(iout)
                lri = maporb(lk+noidx(nsp))
                lrj = maporb(ll+noidx(nsq))
                lrk = maporb(li+noidx(nsr))
                lrl = maporb(lj+noidx(nss))
                !write(u6,'(9(1x,i2),1x,f18.9)') li,lj,lk,ll,lri,lrj,lrk,lrl,iout,val
                call intrw_mol(lri,lrj,lrk,lrl,val)
              end do
            end do
          end do
        end do

      end do
    end do
  end do
end do
call mma_deallocate(buff)

return

2000 format(/7x,'symmetry',6x,' orbitals',8x,'integrals')
2100 format(7x,4i2,1x,4i4,2x,3x,i9)

end subroutine readtwoeint

subroutine readtraonehead(nft,ecor,idisk)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
#include "Molcas.fh"
integer(kind=iwp), intent(in) :: nft
real(kind=wp), intent(out) :: ecor
integer(kind=iwp), intent(out) :: idisk
integer(kind=iwp) :: idum(1), lenrd, nbas(8), ncone(64), ndel(8), nfro(8), norb(8)
real(kind=wp) :: dum(1)
character(len=LenIn8), allocatable :: bsbl(:)

idisk = 0
call idafile(nft,2,ncone,64,idisk)
call ddafile(nft,2,dum,1,idisk)
ecor = dum(1)
call idafile(nft,2,idum,1,idisk)
!nsym = idum(1)
call idafile(nft,2,nbas,8,idisk)
call idafile(nft,2,norb,8,idisk)
call idafile(nft,2,nfro,8,idisk)
call idafile(nft,2,ndel,8,idisk)
call mma_allocate(bsbl,maxbfn,label='bsbl')
lenrd = LenIn8*maxbfn
call cdafile(nft,2,bsbl,lenrd,idisk)
call mma_deallocate(bsbl)

!#ifdef debug
!write(u6,'(a4,1x,8(2x,i8))') 'ncon',ncone(1:8)
!write(u6,*) 'idisk : ', idisk
!write(u6,'(a4,1x,f18.9)') 'ecor',ecor
!write(u6,'(a4,1x,i8)') 'nsym',nsym
!write(u6,'(a4,1x,8(2x,i8))') 'nbas',nbas(1:8)
!write(u6,'(a4,1x,8(2x,i8))') 'norb',norb(1:8)
!write(u6,'(a4,1x,8(2x,i8))') 'nfro',nfro(1:8)
!write(u6,'(a4,1x,8(2x,i8))') 'ndel',ndel(1:8)
!#endif

return

end subroutine readtraonehead

#endif

!subroutine ddatard(nft,xbuff,lend,idisk)
!
!use Definitions, only: wp, iwp
!
!implicit none
!integer(kind=iwp), intent(in) :: nft, lend
!integer(kind=iwp), intent(inout) :: idisk
!real(kind=wp), intent(out) :: xbuff(lend)
!
!lenrd = lend*lendbl
!call idatard(nft,xbuff,lenrd,idisk)
!
!return
!
!end subroutine ddatard

!subroutine idatard(nft,ibuf,lbuf,idisk)
!! write by suo bing, read molcas file
!! imul = 0, the file readed is not a multifile
!! imul = 1, the file readed is a multifile
!
!use Definitions, only: iwp, u6
!
!implicit none
!integer(kind=iwp), intent(in) :: nft, lbuf
!integer(kind=iwp), intent(inout) :: idisk
!character, intent(out) :: ibuf(lbuf)
!integer(kind=iwp) :: imo, ioff, noff, nr
!integer(kind=iwp), parameter :: min_block_length = 512
!
!noff = idisk*min_block_length
!call clseek(nft,noff,nr)
!if (nr /= noff) then
!  write(u6,*) 'error seek file : ',idisk
!  call abend()
!end if
!call cread_molcas(nft,ibuf,lbuf,nr)
!if(nr /= lbuf) then
!  write(u6,*) 'error read file ',nr
!  write(u6,*) 'nft : ',nft
!  write(u6,*) 'lbuf  : ',lbuf
!  write(u6,*) 'idisk : ',idisk
!  call abend()
!else
!  imo = mod(lbuf,min_block_length)
!  if(imo == 0) then
!    ioff = lbuf/min_block_length
!  else
!    ioff = 1+lbuf/min_block_length
!  end if
!  idisk = idisk+ioff
!  !write(u6,*) 'ioff ',ioff
!end if
!
!return
!
!end subroutine idatard

subroutine intrw_mol(ik,jk,kk,lk,val)

use gugaci_global, only: vdint, vector1, voint

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ik, jk, kk, lk
real(kind=wp), intent(in) :: val
integer(kind=iwp) :: i, ind(4), j, k, l, list, lri, lrj, lrk, lrl, lrn, nt
integer(kind=iwp), external :: list3_all, list4_all

list = 0
i = ik
j = jk
k = kk
l = lk
lri = min(i,j)
lrj = max(i,j)
lrk = min(k,l)
lrl = max(k,l)
if (lri > lrk) then
  lrn = lrk
  lrk = lri
  lri = lrn
  lrn = lrl
  lrl = lrj
  lrj = lrn
end if

ind(1) = lri
ind(2) = lrj
ind(3) = lrk
ind(4) = lrl
do i=1,4
  do j=i+1,4
    if (ind(i) > ind(j)) then
      nt = ind(i)
      ind(i) = ind(j)
      ind(j) = nt
    end if
  end do
end do

!write(u6,*) ind(1),ind(2),ind(3),ind(4)

if (ind(1) == ind(4)) then   !(iiii)
  vdint(lri,lri) = val

else if ((lri == lrj) .and. (lrk == lrl)) then   !(iikk)
  vdint(lrk,lri) = val

else if ((lri == lrk) .and. (lrj == lrl)) then   !(ijij)
  !if ((lri == 1) .and. (lrj == 2)) call abend()
  voint(lrj,lri) = val

else if ((lri /= lrl) .and. (lrj == lrk)) then !(ijjl)
  list = list3_all(lri,lrl,lrj)
  vector1(list) = val
  if (lrj == lrl) then
    vector1(list+1) = val
  end if
  if (lri == lrj) then
    vector1(list+1) = val
  end if

else if ((lri /= lrj) .and. (lrk == lrl)) then  !(ijkk)
  list = list3_all(lri,lrj,lrk)
  vector1(list+1) = val
  if (lrj == lrk) then
    vector1(list) = val
  end if

else if ((lri == lrj) .and. (lrk /= lrl)) then  !(iikl)
  list = list3_all(lrk,lrl,lri)
  vector1(list+1) = val

else if ((lri == lrk) .and. (lrj /= lrl)) then      !(ijil)
  if (lrj < lrl) then
    list = list3_all(lrj,lrl,lri)
  else
    list = list3_all(lrl,lrj,lri)
  end if
  vector1(list) = val

else if ((lri /= lrk) .and. (lrj == lrl)) then   !(ijkj)
  list = list3_all(lri,lrk,lrj)
  vector1(list) = val
  if (lrj == lrk) then
    vector1(list+1) = val
  end if

else if ((lri /= lrj) .and. (lri /= lrk) .and. (lri /= lrl) .and. (lrj /= lrk) .and. (lrj /= lrl) .and. (lrk /= lrl)) then
  list = list4_all(ind(1),ind(2),ind(3),ind(4))       !(ijkl)
  if ((lrj > lrk) .and. (lrj > lrl)) then
    vector1(list+2) = val
  end if
  if ((lrj > lrk) .and. (lrj < lrl)) then
    vector1(list) = val
  end if
  if ((lrj < lrk) .and. (lrj < lrl)) then
    vector1(list+1) = val
  end if

end if

!write(u6,*) 'list=',list

return

end subroutine intrw_mol

subroutine int_index(numb)

use gugaci_global, only: loij_all, loijk_all, lsm, ncibl_all, ngw2, ngw3, norb_all, norb_number
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: numb
integer(kind=iwp) :: i, j, la, lb, lc, ld, lra, lrb, lrc, lrd, lri, lrj, ms, msa, msb, msc, mscd, msd, msob(8), nij, njkl

msob = 0
do la=norb_all,1,-1
  lra = norb_all-la+1
  ms = lsm(lra)
  msob(ms) = msob(ms)+1
  ncibl_all(la) = msob(ms)
end do

!write(u6,*)
!=======================================================================
!write(u6,'(1x,14i3)') (ncibl(i),i=1,14)

numb = 1
do i=1,norb_all-1
  do j=i+1,norb_all
    lri = norb_number(i)
    lrj = norb_number(j)
    if (lsm(lri) /= lsm(lrj)) cycle
    !write(u6,*) 'lri=',lri,'lrj=',lrj
    !write(u6,*) 'lsm(lri)',lsm(lri),'lsm(lrj)',lsm(lrj)

    nij = i+ngw2(j)
    !write(u6,*) 'i,j,mij,nij   ',i,j,mij,nij
    loij_all(nij) = numb
    !do k=1,norb_all
    !  vint_ci(numb) = vfutei(j,k,i,k)
    !  vint_ci(numb+1) = vfutei(j,i,k,k)
    !  write(10,'(2x,4i6,i8,3f16.8)') i,j,k,k, numb,vint_ci(numb),vint_ci(numb+1)
    numb = numb+2*norb_all
    !end do

  end do
end do
!write(u6,*) 'number 3 index',numb
!call abend()
!=======================================================================
! la<lb<lc<ld
do ld=1,norb_all-3
  do lc=ld+1,norb_all-2
    lrd = norb_all-ld+1
    lrc = norb_all-lc+1
    msd = lsm(lrd)
    msc = lsm(lrc)
    mscd = Mul(msd,msc)
    do lb=lc+1,norb_all-1
      lrb = norb_all-lb+1
      msb = lsm(lrb)
      msa = Mul(mscd,msb)

      njkl = ld+ngw2(lc)+ngw3(lb)
      loijk_all(njkl) = numb

      !nolra = 0
      do la=norb_all,lb+1,-1
        lra = norb_all-la+1
        if (lsm(lra) /= msa) cycle
        !nolra = nolra+1
        !list = loijk_all(njkl)+3*(nolra-1)
        !write(u6,*) 'ld,lc,lb,la',ld,lc,lb,la,list

        !write(u6,'(2x,4i3,2i7)') la,lb,lc,ld,list,numb

        !vint_ci(numb) = vfutei(la,lc,lb,ld)        !tmp stop
        !vint_ci(numb+1) = vfutei(la,lb,lc,ld)
        !vint_ci(numb+2) = vfutei(la,ld,lc,lb)
        !write(10,'(2x,4i6,i8,3f16.8)') la,lb,lc,ld, numb,vint_ci(numb),vint_ci(numb+1),vint_ci(numb+2)
        numb = numb+3
      end do
    end do
  end do
end do
!call abend()

return

end subroutine int_index

function vfutei(ix,jx,kx,lx)

use gugaci_global, only: vector1
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: vfutei
integer(kind=iwp), intent(in) :: ix, jx, kx, lx
integer(kind=iwp) :: i, ind(4), j, list, lri, lrj, lrk, lrl, lrn, nt
real(kind=wp) :: val
integer(kind=iwp), external :: list3_all, list4_all

lri = min(ix,jx)
lrj = max(ix,jx)
lrk = min(kx,lx)
lrl = max(kx,lx)
if (lri > lrk) then
  lrn = lrk
  lrk = lri
  lri = lrn
  lrn = lrl
  lrl = lrj
  lrj = lrn
end if
val = Zero
!write(u6,*) ind(1),ind(2),ind(3),ind(4)
if ((lri /= lrl) .and. (lrj == lrk)) then !(ijjl)
  list = list3_all(lri,lrl,lrj)
  val = vector1(list)

else if ((lri /= lrj) .and. (lrk == lrl)) then  !(ijkk)
  list = list3_all(lri,lrj,lrk)
  val = vector1(list+1)

else if ((lri == lrj) .and. (lrk /= lrl)) then  !(iikl)
  list = list3_all(lrk,lrl,lri)
  val = vector1(list+1)

else if ((lri == lrk) .and. (lrj /= lrl)) then  !(ijil)
  if (lrj < lrl) then
    list = list3_all(lrj,lrl,lri)
  else
    list = list3_all(lrl,lrj,lri)
  end if
  val = vector1(list)

else if ((lri /= lrk) .and. (lrj == lrl)) then  !(ijkj)
  list = list3_all(lri,lrk,lrj)
  val = vector1(list)

else if ((lri /= lrj) .and. (lri /= lrk) .and. (lri /= lrl) .and. (lrj /= lrk) .and. (lrj /= lrl) .and. (lrk /= lrl)) then
  ind(1) = lri
  ind(2) = lrj
  ind(3) = lrk
  ind(4) = lrl
  do i=1,4
    do j=i+1,4
      if (ind(i) > ind(j)) then
        nt = ind(i)
        ind(i) = ind(j)
        ind(j) = nt
      end if
    end do
  end do

  list = list4_all(ind(1),ind(2),ind(3),ind(4))  !(ijkl)
  if ((lrj > lrk) .and. (lrj > lrl)) then
    val = vector1(list+2)
  end if
  if ((lrj > lrk) .and. (lrj < lrl)) then
    val = vector1(list)
  end if
  if ((lrj < lrk) .and. (lrj < lrl)) then
    val = vector1(list+1)
  end if
end if

vfutei = val

end function vfutei

function list3_all(i,j,k)

use gugaci_global, only: loij_all, ngw2
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: list3_all
integer(kind=iwp), intent(in) :: i, j, k
integer(kind=iwp) :: nij

nij = i+ngw2(j)
list3_all = loij_all(nij)+2*(k-1)

return

end function list3_all

function list4_all(ld,lc,lb,la)

use gugaci_global, only: loijk_all, ncibl_all, ngw2, ngw3
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: list4_all
integer(kind=iwp), intent(in) :: ld, lc, lb, la
integer(kind=iwp) :: lra, njkl

!write(u6,*) 'ld,lc,lb,la',ld,lc,lb,la
lra = ncibl_all(la)
njkl = ld+ngw2(lc)+ngw3(lb)
list4_all = loijk_all(njkl)+3*(lra-1)

return

end function list4_all
