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

subroutine t3loopa(oeh,oep,t1a,t1b,nga,ngb,ngc,vblock,energ,isp,LU,ifvo,scored,enx)
!mp subroutine t3loopa(oeh,oep,t1a,t1b,g,nga,ngb,ngc,vblock,energ,
! implemented integer offsets, PV, 16 may 2004.

use ChT3_global, only: nblock, NNOAB, NOAB, NUAB
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: oeh(*), oep(*)
real(kind=wp), intent(inout) :: t1a(*), t1b(*), energ(*)
integer(kind=iwp), intent(in) :: nga, ngb, ngc, vblock, isp, LU(*)
logical(kind=iwp), intent(in) :: ifvo
logical(kind=iwp), intent(out) :: scored
real(kind=wp), intent(out) :: enx
integer(kind=iwp) :: adim, aset, bdim, bset, cdim, cset, iasblock(3), IUHF, n, nug
real(kind=wp), allocatable :: ka(:), la(:), lb(:), lc(:), t3a(:), t3b(:), voa(:), vob(:), voc(:)
real(kind=wp), allocatable, save :: kb(:), kc(:)

N = noab(isp)+nuab(isp)
enx = Zero
scored = .true.
!mp !!!if (.not. lastcall) then
!mp write(u6,*) 'NOAB,NNOAB,NUAB,NNUAB,ICH'
!mp write(u6,*) NOAB,NNOAB,NUAB,NNUAB,ICH
!mp !!!if (energ(isp) == Zero) then
! this is a first entry - initialization (makes no harm if repeated)
nug = nuab(isp)/vblock
if ((nug*vblock) < nuab(isp)) nug = nug+1
!mp write(u6,*) 'first,nug,vblock',nug,vblock,iopt(1)
IUHF = isp
!!if (IOPT(1) == 0) IUHF = 3
iasblock(1) = vblock*vblock*N/nblock
if (iasblock(1)*nblock < vblock*vblock*N) iasblock(1) = iasblock(1)+1
iasblock(2) = nnoab(iuhf)*vblock*N/nblock
if (iasblock(2)*nblock < nnoab(iuhf)*vblock*N) iasblock(2) = iasblock(2)+1
iasblock(3) = nnoab(iuhf)*vblock*vblock/nblock
if (iasblock(3)*nblock < nnoab(iuhf)*vblock*vblock) iasblock(3) = iasblock(3)+1
!mp call w_rescope(G,'G3loopa')
!mp call w_free(g,0,'G3loopa')
! allocations
!mp call w_alloc(ka,noab(isp)*vblock*vblock*n,'kaT3loopa')
call mma_allocate(ka,noab(isp)*vblock*vblock*n,label='loopa_ka')
if (nug /= 1) then
  !mp call w_alloc(kb,noab(isp)*vblock*vblock*n,'kbT3loopa')
  call mma_allocate(kb,noab(isp)*vblock*vblock*n,label='loopa_kb')
  !mp call w_alloc(kc,noab(isp)*vblock*vblock*n,'kcT3loopa')
  call mma_allocate(kc,noab(isp)*vblock*vblock*n,label='loopa_kc')
end if
!mp call w_alloc(la,nnoab(IUHF)*vblock*n,'laT3loopa')
call mma_allocate(la,nnoab(IUHF)*vblock*n,label='loopa_la')
!mp call w_alloc(lb,nnoab(IUHF)*vblock*n,'lbT3loopa')
call mma_allocate(lb,nnoab(IUHF)*vblock*n,label='loopa_lb')
!mp call w_alloc(lc,nnoab(IUHF)*vblock*n,'lcT3loopa')
call mma_allocate(lc,nnoab(IUHF)*vblock*n,label='loopa_lc')
!mp call w_alloc(t3a,vblock*vblock*vblock,'t3aT3loopa')
call mma_allocate(t3a,vblock*vblock*vblock,label='loopa_t3a')
!mp call w_alloc(t3b,vblock*vblock*vblock,'t3bT3loopa')
call mma_allocate(t3b,vblock*vblock*vblock,label='loopa_t3b')
!mp call w_alloc(voa,vblock*vblock*nnoab(IUHF),'voaT3loopa')
call mma_allocate(voa,vblock*vblock*nnoab(IUHF),label='loopa_voa')
!mp call w_alloc(vob,vblock*vblock*nnoab(IUHF),'vobT3loopa')
call mma_allocate(vob,vblock*vblock*nnoab(IUHF),label='loopa_vob')
!mp call w_alloc(voc,vblock*vblock*nnoab(IUHF),'vocT3loopa')
call mma_allocate(voc,vblock*vblock*nnoab(IUHF),label='loopa_voc')
! this is necessary
! prefactors currently a formal allocation
!mp ? call w_alloc(mi,1,'miT3loopa')
!mp ? call w_alloc(mij,1,'T3loopa')

!mp !!!end if     ! energ - initialization
aset = (nga-1)*vblock
adim = min(vblock,nuab(isp)-aset)
bset = (ngb-1)*vblock
bdim = min(vblock,nuab(isp)-bset)
cset = (ngc-1)*vblock
cdim = min(vblock,nuab(isp)-cset)

! case1 nga=ngb=ngc
! if memory is available loops over i,j,k, in subloops can be grouped !!!
if (nga == ngc) then

  !mp call t3_bt_aaa(nug,g(ka),g(ka),g(ka),g(la),g(mi),g(mij),adim,N,noab(isp),nuab(isp),nnoab(iuhf),lu,iasblock,nga,oeh, &
  !mp                oep(aset+1),enx,g(voa),t1a(noab(isp)*aset+1),t1b(noab(isp)*aset+1),g(t3a),g(t3b),ifvo)
  call t3_bt_aaa(nug,ka,la,adim,N,noab(isp),nnoab(iuhf),lu,iasblock,nga,oeh,oep(aset+1),enx,voa,t1a(noab(isp)*aset+1), &
                 t1b(noab(isp)*aset+1),t3a,t3b,ifvo)

else if (nga == ngb) then
  !mp call t3_bt_aac(nug,g(ka),g(kb),g(kc),g(la),g(lc),g(mi),g(mij),adim,cdim,N,noab(isp),nuab(isp),nnoab(iuhf),lu,iasblock,nga, &
  !mp                ngc,oeh,oep(aset+1),oep(cset+1),enx,g(voa),g(voc),t1a(noab(isp)*aset+1),t1b(noab(isp)*aset+1), &
  !mp                t1a(noab(isp)*cset+1),t1b(noab(isp)*cset+1),g(t3a),g(t3b),ifvo)
  call t3_bt_aac(nug,ka,kc,la,lc,adim,cdim,N,noab(isp),nnoab(iuhf),lu,iasblock,nga,ngc,oeh,oep(aset+1),oep(cset+1),enx,voa,voc, &
                 t1a(noab(isp)*aset+1),t1b(noab(isp)*aset+1),t1a(noab(isp)*cset+1),t1b(noab(isp)*cset+1),t3a,t3b,ifvo)
else if (ngb == ngc) then
  !mp call t3_bt_acc(nug,g(ka),g(kb),g(kc),g(la),g(lc),g(mi),g(mij),adim,cdim,N,noab(isp),nuab(isp),nnoab(iuhf),lu,iasblock,nga, &
  !mp                ngc,oeh,oep(aset+1),oep(cset+1),enx,g(voa),g(voc),t1a(noab(isp)*aset+1),t1b(noab(isp)*aset+1), &
  !mp                t1a(noab(isp)*cset+1),t1b(noab(isp)*cset+1),g(t3a),g(t3b),ifvo)
  call t3_bt_acc(nug,ka,kc,la,lc,adim,cdim,N,noab(isp),nnoab(iuhf),lu,iasblock,nga,ngc,oeh,oep(aset+1),oep(cset+1),enx,voa,voc, &
                 t1a(noab(isp)*aset+1),t1b(noab(isp)*aset+1),t1a(noab(isp)*cset+1),t1b(noab(isp)*cset+1),t3a,t3b,ifvo)
else
  !mp call t3_bt_abc(nug,g(ka),g(kb),g(kc),g(la),g(lb),g(lc),g(mi),g(mij),adim,bdim,cdim,N,noab(isp),nuab(isp),nnoab(iuhf),lu, &
  !mp                iasblock,nga,ngb,ngc,oeh,oep(aset+1),oep(bset+1),oep(cset+1),enx,g(voa),g(vob),g(voc),t1a(noab(isp)*aset+1), &
  !mp                t1b(noab(isp)*aset+1),t1a(noab(isp)*bset+1),t1b(noab(isp)*bset+1),t1a(noab(isp)*cset+1), &
  !mp                t1b(noab(isp)*cset+1),g(t3a),g(t3b),ifvo)
  call t3_bt_abc(nug,ka,kb,kc,la,lb,lc,adim,bdim,cdim,N,noab(isp),nnoab(iuhf),lu,iasblock,nga,ngb,ngc,oeh,oep(aset+1),oep(bset+1), &
                 oep(cset+1),enx,voa,vob,voc,t1a(noab(isp)*aset+1),t1b(noab(isp)*aset+1),t1a(noab(isp)*bset+1), &
                 t1b(noab(isp)*bset+1),t1a(noab(isp)*cset+1),t1b(noab(isp)*cset+1),t3a,t3b,ifvo)
end if   ! cases
energ(isp) = energ(isp)+enx
!mp write(u6,'(A,3(i3,1x),f21.19)') 'nga, ngb, ngc, inc = ',nga,ngb,ngc,enx
!mp
!mp !!!end if  ! lastcall
!mp write(u6,*)
!mp write(u6,*) 'deallocating arrays in t3loopa'
!mp write(u6,*)
call mma_deallocate(voc)
call mma_deallocate(vob)
call mma_deallocate(voa)
call mma_deallocate(t3a)
call mma_deallocate(t3b)
call mma_deallocate(lc)
call mma_deallocate(lb)
call mma_deallocate(la)
if (nug /= 1) then
  call mma_deallocate(kc)
  call mma_deallocate(kb)
end if
call mma_deallocate(ka)

return

end subroutine t3loopa
