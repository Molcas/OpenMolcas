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

subroutine t3loopb(oeh,oep,t1a,t1b,nga,ngb,ngc,vblock,energ,isp,LU,ifvo,scored,enx)
!mp subroutine t3loopb(oeh,oep,t1a,t1b,g,nga,ngb,ngc,vblock,energ,
! implemented integer offsets, PV, 16 may 2004.

use ChT3_global, only: nblock, NNOAB, NNUAB, NOAB, NUAB
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: oeh(*), oep(*)
real(kind=wp), intent(inout) :: t1a(*), t1b(*), energ(*)
integer(kind=iwp), intent(in) :: nga, ngb, ngc, vblock, isp, lu(6)
logical(kind=iwp), intent(in) :: ifvo
logical(kind=iwp), intent(out) :: scored
real(kind=wp), intent(out) :: enx
integer(kind=iwp) :: adim, aset, bdim, bset, cdim, cset, en_offset_ah, en_offset_ap, en_offset_bh, en_offset_bp, iasblock(5), n, &
                     nuga, nugc, t1_offset_a, t1_offset_b
real(kind=wp), allocatable :: kab(:), kac(:), kbc(:), kc(:), kca(:), kcb(:), la(:), lb(:), lxa(:), lxb(:), lxc(:), mi(:), mij(:), &
                              t3a(:), t3b(:), vab(:), vac(:), vbc(:)

en_offset_ah = (isp-1)*noab(1)
en_offset_bh = (2-isp)*noab(1)
en_offset_ap = (isp-1)*nuab(1)
en_offset_bp = (2-isp)*nuab(1)
t1_offset_a = (isp-1)*noab(1)*nuab(1)
t1_offset_b = (2-isp)*noab(1)*nuab(1)
N = noab(isp)+nuab(isp)

!mp
scored = .true.
!mp !!!if (.not. lastcall) then
!mp
enx = Zero
!mpn write(u6,*) 'Bef isp, energ(isp), enx = ',isp,energ(isp),enx
!mp !!!if (energ(isp) == Zero) then
! this is a first entry - initialization (makes no harm if repeated)
nuga = nuab(isp)/vblock
if ((nuga*vblock) < nuab(isp)) nuga = nuga+1
!!write(u6,*) 'first,nuga,vblock',nuga,vblock
nugc = nuab(3-isp)/vblock
if ((nugc*vblock) < nuab(3-isp)) nugc = nugc+1
!!write(u6,*) 'first,nugc,vblock',nugc,vblock

iasblock(1) = vblock*vblock*N/nblock
if (iasblock(1)*nblock < vblock*vblock*N) iasblock(1) = iasblock(1)+1

iasblock(2) = nnoab(isp)*vblock*N/nblock
if (iasblock(2)*nblock < nnoab(isp)*vblock*N) iasblock(2) = iasblock(2)+1

iasblock(3) = nnoab(3)*vblock*N/nblock
if (iasblock(3)*nblock < nnoab(3)*vblock*N) iasblock(3) = iasblock(3)+1

iasblock(4) = nnoab(isp)*vblock*vblock/nblock
if (iasblock(4)*nblock < nnoab(isp)*vblock*vblock) iasblock(4) = iasblock(4)+1

iasblock(5) = nnoab(3)*vblock*vblock/nblock
if (iasblock(5)*nblock < nnoab(3)*vblock*vblock) iasblock(5) = iasblock(5)+1

!mp call w_rescope(G,'GT3loopb')
!mp call w_free(g,0,'GT3loopb')
! allocations
if (nuga /= 1) then
  !mp call w_alloc(kab,noab(isp)*vblock*vblock*n,'kaT3loopb')
  call mma_allocate(kab,noab(isp)*vblock*vblock*n,label='loopb_kab')
  !mp call w_alloc(kcb,noab(isp)*vblock*vblock*n,'kbT3loopb')
  call mma_allocate(kcb,noab(isp)*vblock*vblock*n,label='loopb_kcb')
  !mp call w_alloc(kbc,vblock*vblock*n,'kcT3loopb')
  call mma_allocate(kbc,vblock*vblock*n,label='loopb_kbc')
else
  !mp call w_alloc(kab,noab(isp)*N*nnuab(isp),'kaT3loopb')
  call mma_allocate(kab,noab(isp)*N*nnuab(isp),label='loopb_kab')
end if
!mp call w_alloc(kac,vblock*vblock*n,'kcT3loopb')
call mma_allocate(kac,vblock*vblock*n,label='loopb_kac')
!mp call w_alloc(kca,noab(isp)*vblock*vblock*n,'kbT3loopb')
call mma_allocate(kca,noab(isp)*vblock*vblock*n,label='loopb_kca')
!mp call w_alloc(kc,vblock*vblock*n,'kcT3loopb')
call mma_allocate(kc,vblock*vblock*n,label='loopb_kc')
!mp call w_alloc(la,nnoab(isp)*vblock*n,'laT3loopb')
call mma_allocate(la,nnoab(isp)*vblock*n,label='loopb_la')
!mp call w_alloc(lxa,nnoab(3)*vblock*n,'lbaT3loopb')
call mma_allocate(lxa,nnoab(3)*vblock*n,label='loopb_lxa')
!mpn write(u6,*) 'check lxa'
!mpn write(u6,*) 'nnoab(1), (2), (3) = ',nnoab(1),nnoab(2),nnoab(3)
!mpn write(u6,*) 'nnuab(1), (2), (3) = ',nnuab(1),nnuab(2),nnuab(3)
if (nuga /= 1) then
  !mp call w_alloc(lb,nnoab(isp)*vblock*n,'lbT3loopb')
  call mma_allocate(lb,nnoab(isp)*vblock*n,label='loopb_lb')
  !mp call w_alloc(lxb,nnoab(3)*vblock*n,'labT3loopb')
  call mma_allocate(lxb,nnoab(3)*vblock*n,label='loopb_lxb')
end if
!mp call w_alloc(lxc,nnoab(3)*vblock*n,'lacT3loopb')
call mma_allocate(lxc,nnoab(3)*vblock*n,label='loopb_lxc')
!mp call w_alloc(t3a,vblock*vblock*vblock,'t3aT3loopb')
call mma_allocate(t3a,vblock*vblock*vblock,label='loopb_t3a')
!mp call w_alloc(t3b,vblock*vblock*vblock,'t3bT3loopb')
call mma_allocate(t3b,vblock*vblock*vblock,label='loopb_t3b')
!mp call w_alloc(vac,vblock*vblock*nnoab(3),'vbcT3loopb')
call mma_allocate(vac,vblock*vblock*nnoab(3),label='loopb_vac')
if (nuga /= 1) then
  !mp call w_alloc(vab,vblock*vblock*nnoab(isp),'vabT3loopb')
  call mma_allocate(vab,vblock*vblock*nnoab(isp),label='loopb_vab')
  !mp call w_alloc(vbc,vblock*vblock*nnoab(3),'vacT3loopb')
  call mma_allocate(vbc,vblock*vblock*nnoab(3),label='loopb_vbc')
else
  !mp call w_alloc(vab,nnoab(isp)*nnuab(isp),'vabT3loopb')
  call mma_allocate(vab,nnoab(isp)*nnuab(isp),label='loopb_vab')
end if
!mp call w_alloc(mi,noab(isp)*vblock**3,'miT3loopb')
call mma_allocate(mi,noab(isp)*vblock**3,label='loopb_mi')
!mp
!mpn write(u6,*) 'mi check = ',noab(isp)*vblock**3
!mpn write(u6,*) 'nnoab(1), (2), (3) = ',nnoab(1),nnoab(2),nnoab(3)
!mpn write(u6,*) 'nnuab(1), (2), (3) = ',nnuab(1),nnuab(2),nnuab(3)
!mp
!mp call w_alloc(mij,N*vblock,'mijT3loopb')
call mma_allocate(mij,N*vblock,label='loopb_mij')
!mp !!!end if     ! energ = 0 - initialization

aset = (nga-1)*vblock
adim = min(vblock,nuab(isp)-aset)
bset = (ngb-1)*vblock
bdim = min(vblock,nuab(isp)-bset)
cset = (ngc-1)*vblock
cdim = min(vblock,nuab(3-isp)-cset)

! case1 nga=ngb=ngc

if (nga == ngb) then

  !mpn write(u6,*) 'ide call t3_bta_aac'
  !mp call t3_bta_aac(nuga,nugc,g(kab),g(kca),g(kac),g(kc),g(la),g(lxa),g(lxc),g(mi),g(mij),adim,cdim,N,noab(isp),nuab(isp), &
  !mp                 noab(3-isp),nuab(3-isp),lu,iasblock,nga,ngc,oeh(en_offset_ah+1),oeh(en_offset_bh+1), &
  !mp                 oep(aset+en_offset_ap+1),oep(cset+en_offset_bp+1),enx,g(vab),g(vac),t1a(noab(isp)*aset+t1_offset_a+1), &
  !mp                 t1b(noab(isp)*aset+t1_offset_a+1),t1a(noab(3-isp)*cset+t1_offset_b+1),t1b(noab(3-isp)*cset+t1_offset_b+1), &
  !mp                 g(t3a),g(t3b),ifvo)
  call t3_bta_aac(nuga,nugc,kab,kca,kac,kc,la,lxa,lxc,mi,mij,adim,cdim,N,noab(isp),noab(3-isp),lu,iasblock,nga,ngc, &
                  oeh(en_offset_ah+1),oeh(en_offset_bh+1),oep(aset+en_offset_ap+1),oep(cset+en_offset_bp+1),enx,vab,vac, &
                  t1a(noab(isp)*aset+t1_offset_a+1),t1b(noab(isp)*aset+t1_offset_a+1),t1a(noab(3-isp)*cset+t1_offset_b+1), &
                  t1b(noab(3-isp)*cset+t1_offset_b+1),t3a,t3b,ifvo)
else
  !mpn write(u6,*) 'ide call t3_bta_abc'
  !mp call t3_bta_abc(nuga,nugc,g(kab),g(kcb),g(kca),g(kac),g(kbc),g(kc),g(la),g(lb),g(lxa),g(lxb),g(lxc),g(mi),g(mij),adim,bdim, &
  !mp                 cdim,N,noab(isp),nuab(isp),noab(3-isp),nuab(3-isp),lu,iasblock,nga,ngb,ngc,oeh(en_offset_ah+1), &
  !mp                 oeh(en_offset_bh+1),oep(aset+en_offset_ap+1),oep(bset+en_offset_ap+1),oep(cset+en_offset_bp+1),enx,g(vab), &
  !mp                 g(vbc),g(vac),t1a(noab(isp)*aset+t1_offset_a+1),t1b(noab(isp)*aset+t1_offset_a+1), &
  !mp                 t1a(noab(isp)*bset+t1_offset_a+1),t1b(noab(isp)*bset+t1_offset_a+1),t1a(noab(3-isp)*cset+t1_offset_b+1), &
  !mp                 t1b(noab(3-isp)*cset+t1_offset_b+1),g(t3a),g(t3b),ifvo)
  call t3_bta_abc(nuga,nugc,kab,kcb,kca,kac,kbc,kc,la,lb,lxa,lxb,lxc,mi,mij,adim,bdim,cdim,N,noab(isp),noab(3-isp),lu,iasblock, &
                  nga,ngb,ngc,oeh(en_offset_ah+1),oeh(en_offset_bh+1),oep(aset+en_offset_ap+1),oep(bset+en_offset_ap+1), &
                  oep(cset+en_offset_bp+1),enx,vab,vbc,vac,t1a(noab(isp)*aset+t1_offset_a+1),t1b(noab(isp)*aset+t1_offset_a+1), &
                  t1a(noab(isp)*bset+t1_offset_a+1),t1b(noab(isp)*bset+t1_offset_a+1),t1a(noab(3-isp)*cset+t1_offset_b+1), &
                  t1b(noab(3-isp)*cset+t1_offset_b+1),t3a,t3b,ifvo)
end if   ! cases
!mpn write(u6,*) 'isp, energ(isp), enx = ',isp,energ(isp),enx
energ(isp) = energ(isp)+enx
!mp !!!write(u6,'(A,3(i5,2x),f21.19)') 'nga, ngb, ngc, inc = ',nga,ngb,ngc,enx
!mp !!!end if  ! lastcall
!mp write(u6,*)
!mp write(u6,*) 'deallocating arrays in t3loob'
!mp write(u6,*)
!mp
call mma_deallocate(mij)
call mma_deallocate(mi)
if (nuga /= 1) then
  call mma_deallocate(vbc)
end if
call mma_deallocate(vab)
call mma_deallocate(vac)
call mma_deallocate(t3a)
call mma_deallocate(t3b)
call mma_deallocate(lxc)
if (nuga /= 1) then
  call mma_deallocate(lxb)
  call mma_deallocate(lb)
end if
call mma_deallocate(lxa)
call mma_deallocate(la)
call mma_deallocate(kc)
call mma_deallocate(kca)
call mma_deallocate(kac)
if (nuga /= 1) then
  call mma_deallocate(kbc)
  call mma_deallocate(kcb)
end if
call mma_deallocate(kab)

return

end subroutine t3loopb
