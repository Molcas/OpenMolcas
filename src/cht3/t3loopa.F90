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

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: oeh(*), oep(*), t1a(*), t1b(*), energ(*), enx
integer(kind=iwp) :: nga, ngb, ngc, vblock, isp, LU(*)
logical(kind=iwp) :: ifvo, scored
#include "WrkSpc.fh"
#include "ndisk.fh"
#include "uhf.fh"
integer(kind=iwp) :: adim, aset, bdim, bset, cdim, cset, iasblock(3), IUHF, ka, la, lb, lc, mi, mij, n, nug, t3a, t3b, voa, vob, voc
integer(kind=iwp), save :: kb, kc

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
!mp write(u6,*) 'first,nug,vblock',nug,vblock,iopt(76)
IUHF = isp
!!if (IOPT(76) == 0) IUHF = 3
iasblock(1) = vblock*vblock*N/nblock
if ((iasblock(1)*nblock) < (vblock*vblock*N)) iasblock(1) = iasblock(1)+1
iasblock(2) = nnoab(iuhf)*vblock*N/nblock
if ((iasblock(2)*nblock) < (nnoab(iuhf)*vblock*N)) iasblock(2) = iasblock(2)+1
iasblock(3) = nnoab(iuhf)*vblock*vblock/nblock
if ((iasblock(3)*nblock) < (nnoab(iuhf)*vblock*vblock)) iasblock(3) = iasblock(3)+1
!mp call w_rescope(G,'G3loopa')
!mp call w_free(g,0,'G3loopa')
! allocations
!mp call w_alloc(ka,noab(isp)*vblock*vblock*n,'kaT3loopa')
call GetMem('loopa_ka','Allo','Real',ka,noab(isp)*vblock*vblock*n)
if (nug /= 1) then
  !mp call w_alloc(kb,noab(isp)*vblock*vblock*n,'kbT3loopa')
  call GetMem('loopa_kb','Allo','Real',kb,noab(isp)*vblock*vblock*n)
  !mp call w_alloc(kc,noab(isp)*vblock*vblock*n,'kcT3loopa')
  call GetMem('loopa_kc','Allo','Real',kc,noab(isp)*vblock*vblock*n)
end if
!mp call w_alloc(la,nnoab(IUHF)*vblock*n,'laT3loopa')
call GetMem('loopa_la','Allo','Real',la,nnoab(IUHF)*vblock*n)
!mp call w_alloc(lb,nnoab(IUHF)*vblock*n,'lbT3loopa')
call GetMem('loopa_lb','Allo','Real',lb,nnoab(IUHF)*vblock*n)
!mp call w_alloc(lc,nnoab(IUHF)*vblock*n,'lcT3loopa')
call GetMem('loopa_lc','Allo','Real',lc,nnoab(IUHF)*vblock*n)
!mp call w_alloc(t3a,vblock*vblock*vblock,'t3aT3loopa')
call GetMem('loopa_t3a','Allo','Real',t3a,vblock*vblock*vblock)
!mp call w_alloc(t3b,vblock*vblock*vblock,'t3bT3loopa')
call GetMem('loopa_t3b','Allo','Real',t3b,vblock*vblock*vblock)
!mp call w_alloc(voa,vblock*vblock*nnoab(IUHF),'voaT3loopa')
call GetMem('loopa_voa','Allo','Real',voa,vblock*vblock*nnoab(IUHF))
!mp call w_alloc(vob,vblock*vblock*nnoab(IUHF),'vobT3loopa')
call GetMem('loopa_vob','Allo','Real',vob,vblock*vblock*nnoab(IUHF))
!mp call w_alloc(voc,vblock*vblock*nnoab(IUHF),'vocT3loopa')
call GetMem('loopa_voc','Allo','Real',voc,vblock*vblock*nnoab(IUHF))
! this is necessary
! prefactors currently a formal allocation
!mp ? call w_alloc(mi,1,'miT3loopa')
!mp ? call w_alloc(mij,1,'T3loopa')
call GetMem('loopa_mi','Allo','Real',mi,1)
call GetMem('loopa_mij','Allo','Real',mij,1)

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
  call t3_bt_aaa(nug,Work(ka),Work(la),adim,N,noab(isp),nnoab(iuhf),lu,iasblock,nga,oeh,oep(aset+1),enx,Work(voa), &
                 t1a(noab(isp)*aset+1),t1b(noab(isp)*aset+1),Work(t3a),Work(t3b),ifvo)

else if (nga == ngb) then
  !mp call t3_bt_aac(nug,g(ka),g(kb),g(kc),g(la),g(lc),g(mi),g(mij),adim,cdim,N,noab(isp),nuab(isp),nnoab(iuhf),lu,iasblock,nga, &
  !mp                ngc,oeh,oep(aset+1),oep(cset+1),enx,g(voa),g(voc),t1a(noab(isp)*aset+1),t1b(noab(isp)*aset+1), &
  !mp                t1a(noab(isp)*cset+1),t1b(noab(isp)*cset+1),g(t3a),g(t3b),ifvo)
  call t3_bt_aac(nug,Work(ka),Work(kc),Work(la),Work(lc),adim,cdim,N,noab(isp),nnoab(iuhf),lu,iasblock,nga,ngc,oeh, &
                 oep(aset+1),oep(cset+1),enx,Work(voa),Work(voc),t1a(noab(isp)*aset+1),t1b(noab(isp)*aset+1), &
                 t1a(noab(isp)*cset+1),t1b(noab(isp)*cset+1),Work(t3a),Work(t3b),ifvo)
else if (ngb == ngc) then
  !mp call t3_bt_acc(nug,g(ka),g(kb),g(kc),g(la),g(lc),g(mi),g(mij),adim,cdim,N,noab(isp),nuab(isp),nnoab(iuhf),lu,iasblock,nga, &
  !mp                ngc,oeh,oep(aset+1),oep(cset+1),enx,g(voa),g(voc),t1a(noab(isp)*aset+1),t1b(noab(isp)*aset+1), &
  !mp                t1a(noab(isp)*cset+1),t1b(noab(isp)*cset+1),g(t3a),g(t3b),ifvo)
  call t3_bt_acc(nug,Work(ka),Work(kc),Work(la),Work(lc),adim,cdim,N,noab(isp),nnoab(iuhf),lu,iasblock,nga,ngc,oeh,oep(aset+1), &
                 oep(cset+1),enx,Work(voa),Work(voc),t1a(noab(isp)*aset+1),t1b(noab(isp)*aset+1),t1a(noab(isp)*cset+1), &
                 t1b(noab(isp)*cset+1),Work(t3a),Work(t3b),ifvo)
else
  !mp call t3_bt_abc(nug,g(ka),g(kb),g(kc),g(la),g(lb),g(lc),g(mi),g(mij),adim,bdim,cdim,N,noab(isp),nuab(isp),nnoab(iuhf),lu, &
  !mp                iasblock,nga,ngb,ngc,oeh,oep(aset+1),oep(bset+1),oep(cset+1),enx,g(voa),g(vob),g(voc),t1a(noab(isp)*aset+1), &
  !mp                t1b(noab(isp)*aset+1),t1a(noab(isp)*bset+1),t1b(noab(isp)*bset+1),t1a(noab(isp)*cset+1), &
  !mp                t1b(noab(isp)*cset+1),g(t3a),g(t3b),ifvo)
  call t3_bt_abc(nug,Work(ka),Work(kb),Work(kc),Work(la),Work(lb),Work(lc),adim,bdim,cdim,N,noab(isp),nnoab(iuhf),lu,iasblock,nga, &
                 ngb,ngc,oeh,oep(aset+1),oep(bset+1),oep(cset+1),enx,Work(voa),Work(vob),Work(voc),t1a(noab(isp)*aset+1), &
                 t1b(noab(isp)*aset+1),t1a(noab(isp)*bset+1),t1b(noab(isp)*bset+1),t1a(noab(isp)*cset+1),t1b(noab(isp)*cset+1), &
                 Work(t3a),Work(t3b),ifvo)
end if   ! cases
energ(isp) = energ(isp)+enx
!mp write(u6,'(A,3(i3,1x),f21.19)') 'nga, ngb, ngc, inc = ',nga,ngb,ngc,enx
!mp
!mp !!!end if  ! lastcall
!mp write(u6,*)
!mp write(u6,*) 'deallocating arrays in t3loopa'
!mp write(u6,*)
call GetMem('loopa_mij','Free','Real',mij,1)
call GetMem('loopa_mi','Free','Real',mi,1)
call GetMem('loopa_voc','Free','Real',voc,vblock*vblock*nnoab(IUHF))
call GetMem('loopa_vob','Free','Real',vob,vblock*vblock*nnoab(IUHF))
call GetMem('loopa_voa','Free','Real',voa,vblock*vblock*nnoab(IUHF))
call GetMem('loopa_t3b','Free','Real',t3b,vblock*vblock*vblock)
call GetMem('loopa_t3a','Free','Real',t3a,vblock*vblock*vblock)
call GetMem('loopa_lc','Free','Real',lc,nnoab(IUHF)*vblock*n)
call GetMem('loopa_lb','Free','Real',lb,nnoab(IUHF)*vblock*n)
call GetMem('loopa_la','Free','Real',la,nnoab(IUHF)*vblock*n)
if (nug /= 1) then
  call GetMem('loopa_kc','Free','Real',kc,noab(isp)*vblock*vblock*n)
  call GetMem('loopa_kb','Free','Real',kb,noab(isp)*vblock*vblock*n)
end if
call GetMem('loopa_ka','Free','Real',ka,noab(isp)*vblock*vblock*n)

return

end subroutine t3loopa
