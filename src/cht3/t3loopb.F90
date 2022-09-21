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

implicit none
#include "ndisk.fh"
#include "WrkSpc.fh"
!mp real*8 g(*), energ(*), oeh(*), oep(*), enx, t1a(*), t1b(*)
real*8 energ(*), oeh(*), oep(*), enx, t1a(*), t1b(*)
logical ifvo, scored
integer isp, vblock, n, lu(6), nga, ngb, ngc, adim, bdim, cdim
integer en_offset_ah, en_offset_bh, en_offset_ap, en_offset_bp
integer t1_offset_a, t1_offset_b
integer nuga, nugc
integer iasblock(5), aset, bset, cset
#include "uhf.fh"
#include "ioind.fh"
integer kab, kca, kcb, kac, kbc, kc, la, lb, lxa, lxb, lxc, t3a, t3b, vab, vbc, vac, mij, mi
!mp save kab, kca, kcb, kac, kbc, kc, la, lb, lxa, lxb, lxc, t3a, t3b, vab, vbc, vac, mij, mi, iasblock, nuga, nugc

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
enx = 0.d0
!mpn write(6,*) 'Bef isp, energ(isp), enx = ',isp,energ(isp),enx
!mp !!!if (energ(isp) == 0.d0) then
! this is a first entry - initialization (makes no harm if repeated)
nuga = nuab(isp)/vblock
if ((nuga*vblock) < nuab(isp)) nuga = nuga+1
!!write(6,*) 'first,nuga,vblock',nuga,vblock
nugc = nuab(3-isp)/vblock
if ((nugc*vblock) < nuab(3-isp)) nugc = nugc+1
!!write(6,*) 'first,nugc,vblock',nugc,vblock

iasblock(1) = vblock*vblock*N/nblock
if ((iasblock(1)*nblock) < (vblock*vblock*N)) iasblock(1) = iasblock(1)+1

iasblock(2) = nnoab(isp)*vblock*N/nblock
if ((iasblock(2)*nblock) < (nnoab(isp)*vblock*N)) iasblock(2) = iasblock(2)+1

iasblock(3) = nnoab(3)*vblock*N/nblock
if ((iasblock(3)*nblock) < (nnoab(3)*vblock*N)) iasblock(3) = iasblock(3)+1

iasblock(4) = nnoab(isp)*vblock*vblock/nblock
if ((iasblock(4)*nblock) < (nnoab(isp)*vblock*vblock)) iasblock(4) = iasblock(4)+1

iasblock(5) = nnoab(3)*vblock*vblock/nblock
if ((iasblock(5)*nblock) < (nnoab(3)*vblock*vblock)) iasblock(5) = iasblock(5)+1

!mp call w_rescope(G,'GT3loopb')
!mp call w_free(g,0,'GT3loopb')
! allocations
if (nuga /= 1) then
  !mp call w_alloc(kab,noab(isp)*vblock*vblock*n,'kaT3loopb')
  call GetMem('loopb_kab','Allo','Real',kab,noab(isp)*vblock*vblock*n)
  !mp call w_alloc(kcb,noab(isp)*vblock*vblock*n,'kbT3loopb')
  call GetMem('loopb_kcb','Allo','Real',kcb,noab(isp)*vblock*vblock*n)
  !mp call w_alloc(kbc,vblock*vblock*n,'kcT3loopb')
  call GetMem('loopb_kbc','Allo','Real',kbc,vblock*vblock*n)
else
  !mp call w_alloc(kab,noab(isp)*N*nnuab(isp),'kaT3loopb')
  call GetMem('loopb_kab','Allo','Real',kab,noab(isp)*N*nnuab(isp))
end if
!mp call w_alloc(kac,vblock*vblock*n,'kcT3loopb')
call GetMem('loopb_kac','Allo','Real',kac,vblock*vblock*n)
!mp call w_alloc(kca,noab(isp)*vblock*vblock*n,'kbT3loopb')
call GetMem('loopb_kca','Allo','Real',kca,noab(isp)*vblock*vblock*n)
!mp call w_alloc(kc,vblock*vblock*n,'kcT3loopb')
call GetMem('loopb_kc','Allo','Real',kc,vblock*vblock*n)
!mp call w_alloc(la,nnoab(isp)*vblock*n,'laT3loopb')
call GetMem('loopb_la','Allo','Real',la,nnoab(isp)*vblock*n)
!mp call w_alloc(lxa,nnoab(3)*vblock*n,'lbaT3loopb')
call GetMem('loopb_lxa','Allo','Real',lxa,nnoab(3)*vblock*n)
!mpn write(6,*) 'check lxa'
!mpn write(6,*) 'nnoab(1), (2), (3) = ',nnoab(1),nnoab(2),nnoab(3)
!mpn write(6,*) 'nnuab(1), (2), (3) = ',nnuab(1),nnuab(2),nnuab(3)
if (nuga /= 1) then
  !mp call w_alloc(lb,nnoab(isp)*vblock*n,'lbT3loopb')
  call GetMem('loopb_lb','Allo','Real',lb,nnoab(isp)*vblock*n)
  !mp call w_alloc(lxb,nnoab(3)*vblock*n,'labT3loopb')
  call GetMem('loopb_lxb','Allo','Real',lxb,nnoab(3)*vblock*n)
end if
!mp call w_alloc(lxc,nnoab(3)*vblock*n,'lacT3loopb')
call GetMem('loopb_lxc','Allo','Real',lxc,nnoab(3)*vblock*n)
!mp call w_alloc(t3a,vblock*vblock*vblock,'t3aT3loopb')
call GetMem('loopb_t3a','Allo','Real',t3a,vblock*vblock*vblock)
!mp call w_alloc(t3b,vblock*vblock*vblock,'t3bT3loopb')
call GetMem('loopb_t3b','Allo','Real',t3b,vblock*vblock*vblock)
!mp call w_alloc(vac,vblock*vblock*nnoab(3),'vbcT3loopb')
call GetMem('loopb_vac','Allo','Real',vac,vblock*vblock*nnoab(3))
if (nuga /= 1) then
  !mp call w_alloc(vab,vblock*vblock*nnoab(isp),'vabT3loopb')
  call GetMem('loopb_vab','Allo','Real',vab,vblock*vblock*nnoab(isp))
  !mp call w_alloc(vbc,vblock*vblock*nnoab(3),'vacT3loopb')
  call GetMem('loopb_vbc','Allo','Real',vbc,vblock*vblock*nnoab(3))
else
  !mp call w_alloc(vab,nnoab(isp)*nnuab(isp),'vabT3loopb')
  call GetMem('loopb_vab','Allo','Real',vab,nnoab(isp)*nnuab(isp))
end if
!mp call w_alloc(mi,noab(isp)*vblock**3,'miT3loopb')
call GetMem('loopb_mi','Allo','Real',mi,noab(isp)*vblock**3)
!mp
!mpn write(6,*) 'mi check = ',noab(isp)*vblock**3
!mpn write(6,*) 'nnoab(1), (2), (3) = ',nnoab(1),nnoab(2),nnoab(3)
!mpn write(6,*) 'nnuab(1), (2), (3) = ',nnuab(1),nnuab(2),nnuab(3)
!mp
!mp call w_alloc(mij,N*vblock,'mijT3loopb')
call GetMem('loopb_mij','Allo','Real',mij,N*vblock)
!mp !!!end if     ! energ = 0 - initialization

aset = (nga-1)*vblock
adim = min(vblock,nuab(isp)-aset)
bset = (ngb-1)*vblock
bdim = min(vblock,nuab(isp)-bset)
cset = (ngc-1)*vblock
cdim = min(vblock,nuab(3-isp)-cset)

! case1 nga=ngb=ngc

if (nga == ngb) then

  !mpn write(6,*) 'ide call t3_bta_aac'
  !mp call t3_bta_aac(nuga,nugc,g(kab),g(kca),g(kac),g(kc),g(la),g(lxa),g(lxc),g(mi),g(mij),adim,cdim,N,noab(isp),nuab(isp), &
  !mp                 noab(3-isp),nuab(3-isp),lu,iasblock,nga,ngc,oeh(en_offset_ah+1),oeh(en_offset_bh+1), &
  !mp                 oep(aset+en_offset_ap+1),oep(cset+en_offset_bp+1),enx,g(vab),g(vac),t1a(noab(isp)*aset+t1_offset_a+1), &
  !mp                 t1b(noab(isp)*aset+t1_offset_a+1),t1a(noab(3-isp)*cset+t1_offset_b+1),t1b(noab(3-isp)*cset+t1_offset_b+1), &
  !mp                 g(t3a),g(t3b),ifvo)
  call t3_bta_aac(nuga,nugc,Work(kab),Work(kca),Work(kac),Work(kc),Work(la),Work(lxa),Work(lxc),Work(mi),Work(mij),adim,cdim,N, &
                  noab(isp),noab(3-isp),lu,iasblock,nga,ngc,oeh(en_offset_ah+1),oeh(en_offset_bh+1),oep(aset+en_offset_ap+1), &
                  oep(cset+en_offset_bp+1),enx,Work(vab),Work(vac),t1a(noab(isp)*aset+t1_offset_a+1), &
                  t1b(noab(isp)*aset+t1_offset_a+1),t1a(noab(3-isp)*cset+t1_offset_b+1),t1b(noab(3-isp)*cset+t1_offset_b+1), &
                  Work(t3a),Work(t3b),ifvo)
else
  !mpn write(6,*) 'ide call t3_bta_abc'
  !mp call t3_bta_abc(nuga,nugc,g(kab),g(kcb),g(kca),g(kac),g(kbc),g(kc),g(la),g(lb),g(lxa),g(lxb),g(lxc),g(mi),g(mij),adim,bdim, &
  !mp                 cdim,N,noab(isp),nuab(isp),noab(3-isp),nuab(3-isp),lu,iasblock,nga,ngb,ngc,oeh(en_offset_ah+1), &
  !mp                 oeh(en_offset_bh+1),oep(aset+en_offset_ap+1),oep(bset+en_offset_ap+1),oep(cset+en_offset_bp+1),enx,g(vab), &
  !mp                 g(vbc),g(vac),t1a(noab(isp)*aset+t1_offset_a+1),t1b(noab(isp)*aset+t1_offset_a+1), &
  !mp                 t1a(noab(isp)*bset+t1_offset_a+1),t1b(noab(isp)*bset+t1_offset_a+1),t1a(noab(3-isp)*cset+t1_offset_b+1), &
  !mp                 t1b(noab(3-isp)*cset+t1_offset_b+1),g(t3a),g(t3b),ifvo)
  call t3_bta_abc(nuga,nugc,Work(kab),Work(kcb),Work(kca),Work(kac),Work(kbc),Work(kc),Work(la),Work(lb),Work(lxa),Work(lxb), &
                  Work(lxc),Work(mi),Work(mij),adim,bdim,cdim,N,noab(isp),noab(3-isp),lu,iasblock,nga,ngb,ngc,oeh(en_offset_ah+1), &
                  oeh(en_offset_bh+1),oep(aset+en_offset_ap+1),oep(bset+en_offset_ap+1),oep(cset+en_offset_bp+1),enx,Work(vab), &
                  Work(vbc),Work(vac),t1a(noab(isp)*aset+t1_offset_a+1),t1b(noab(isp)*aset+t1_offset_a+1), &
                  t1a(noab(isp)*bset+t1_offset_a+1),t1b(noab(isp)*bset+t1_offset_a+1),t1a(noab(3-isp)*cset+t1_offset_b+1), &
                  t1b(noab(3-isp)*cset+t1_offset_b+1),Work(t3a),Work(t3b),ifvo)
end if   ! cases
!mpn write(6,*) 'isp, energ(isp), enx = ',isp,energ(isp),enx
energ(isp) = energ(isp)+enx
!mp !!!write(6,'(A,3(i5,2x),f21.19)') 'nga, ngb, ngc, inc = ',nga,ngb,ngc,enx
!mp !!!end if  ! lastcall
!mp write(6,*)
!mp write(6,*) 'deallocating arrays in t3loob'
!mp write(6,*)
!mp
call GetMem('loopb_mij','Free','Real',mij,N*vblock)
call GetMem('loopb_mi','Free','Real',mi,noab(isp)*vblock**3)
if (nuga /= 1) then
  call GetMem('loopb_vbc','Free','Real',vbc,vblock*vblock*nnoab(3))
  call GetMem('loopb_vab','Free','Real',vab,vblock*vblock*nnoab(isp))
else
  call GetMem('loopb_vab','Free','Real',vab,nnoab(isp)*nnuab(isp))
end if
call GetMem('loopb_vac','Free','Real',vac,vblock*vblock*nnoab(3))
call GetMem('loopb_t3b','Free','Real',t3b,vblock*vblock*vblock)
call GetMem('loopb_t3a','Free','Real',t3a,vblock*vblock*vblock)
call GetMem('loopb_lxc','Free','Real',lxc,nnoab(3)*vblock*n)
if (nuga /= 1) then
  call GetMem('loopb_lxb','Free','Real',lxb,nnoab(3)*vblock*n)
  call GetMem('loopb_lb','Free','Real',lb,nnoab(isp)*vblock*n)
end if
call GetMem('loopb_lxa','Free','Real',lxa,nnoab(3)*vblock*n)
call GetMem('loopb_la','Free','Real',la,nnoab(isp)*vblock*n)
call GetMem('loopb_kc','Free','Real',kc,vblock*vblock*n)
call GetMem('loopb_kca','Free','Real',kca,noab(isp)*vblock*vblock*n)
call GetMem('loopb_kac','Free','Real',kac,vblock*vblock*n)
if (nuga /= 1) then
  call GetMem('loopb_kbc','Free','Real',kbc,vblock*vblock*n)
  call GetMem('loopb_kcb','Free','Real',kcb,noab(isp)*vblock*vblock*n)
  call GetMem('loopb_kab','Free','Real',kab,noab(isp)*vblock*vblock*n)
else
  call GetMem('loopb_kab','Free','Real',kab,noab(isp)*N*nnuab(isp))
end if

return

end subroutine t3loopb
