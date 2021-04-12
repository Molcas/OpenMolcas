************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE t3loopa(oeh,oep,t1a,t1b,nga,ngb,ngc,vblock,energ,
     $ isp,LU,ifvo,lastcall,scored,jjj,enx)
cmp      SUBROUTINE t3loopa(oeh,oep,t1a,t1b,g,nga,ngb,ngc,vblock,energ,
C  implemented integer offsets, PV, 16 may 2004.
      IMPLICIT NONE
#include "ndisk.fh"
#include "WrkSpc.fh"
cmp      real*8 g(*),energ(*),oeh(*),oep(*),enx,t1a(*),t1b(*)
      real*8 energ(*),oeh(*),oep(*),enx,t1a(*),t1b(*)
      integer nug
      integer isp,vblock,n,lu(*),nga,ngb,ngc,adim,bdim,cdim
      INTEGER iasblock(3),aset,bset,cset
      logical ifvo,lastcall,scored
      INTEGER IUHF
cmp
        integer jjj
cmp
#include "uhf.fh"
#include "ioind.fh"
      integer ka,kb,kc,la,lb,lc,t3a,t3b,voa,vob,voc,mi,mij
      SAVE      ka,kb,kc,la,lb,lc,t3a,t3b,voa,vob,voc,mi,mij,
     $     iasblock,iuhf,nug
C
      N=noab(isp)+nuab(isp)
      enx=0.d0
        scored=.true.
cmp!!!  if (lastcall) goto 321
cmp     write (6,*) 'NOAB,NNOAB,NUAB,NNUAB,ICH'
cmp     write (6,*) NOAB,NNOAB,NUAB,NNUAB,ICH
cmp!!!      if(energ(isp).eq.0.d0)then
C this is a first entry - initialization (makes no harm if reapeated)
      nug=nuab(isp)/vblock
      if((nug*vblock).lt.nuab(isp))nug=nug+1
cmp      write(6,*)'first,nug,vblock',nug,vblock,iopt(76)
      IUHF=isp
      !!IF(IOPT(76).eq.0)IUHF=3
      iasblock(1)=vblock*vblock*N/nblock
      if((iasblock(1)*nblock).lt.(vblock*vblock*N))
     $                                    iasblock(1)=iasblock(1)+1
      iasblock(2)=nnoab(iuhf)*vblock*N/nblock
      if((iasblock(2)*nblock).lt.(nnoab(iuhf)*vblock*N))
     $iasblock(2)=iasblock(2)+1
      iasblock(3)=nnoab(iuhf)*vblock*vblock/nblock
      if((iasblock(3)*nblock).lt.(nnoab(iuhf)*vblock*vblock))
     $iasblock(3)=iasblock(3)+1
cmp      call w_rescope(G,'G3loopa')
cmp      call w_free(g,0,'G3loopa')
c  allocations
cmp      call w_alloc(ka,noab(isp)*vblock*vblock*n,'kaT3loopa')
      call GetMem('loopa_ka','Allo','Real',ka,noab(isp)*vblock*vblock*n)
      if(nug.ne.1)then
cmp      call w_alloc(kb,noab(isp)*vblock*vblock*n,'kbT3loopa')
      call GetMem('loopa_kb','Allo','Real',kb,noab(isp)*vblock*vblock*n)
cmp      call w_alloc(kc,noab(isp)*vblock*vblock*n,'kcT3loopa')
      call GetMem('loopa_kc','Allo','Real',kc,noab(isp)*vblock*vblock*n)
      endif
cmp      call w_alloc(la,nnoab(IUHF)*vblock*n,'laT3loopa')
      call GetMem('loopa_la','Allo','Real',la,nnoab(IUHF)*vblock*n)
cmp      call w_alloc(lb,nnoab(IUHF)*vblock*n,'lbT3loopa')
      call GetMem('loopa_lb','Allo','Real',lb,nnoab(IUHF)*vblock*n)
cmp      call w_alloc(lc,nnoab(IUHF)*vblock*n,'lcT3loopa')
      call GetMem('loopa_lc','Allo','Real',lc,nnoab(IUHF)*vblock*n)
cmp      call w_alloc(t3a,vblock*vblock*vblock,'t3aT3loopa')
      call GetMem('loopa_t3a','Allo','Real',t3a,vblock*vblock*vblock)
cmp      call w_alloc(t3b,vblock*vblock*vblock,'t3bT3loopa')
      call GetMem('loopa_t3b','Allo','Real',t3b,vblock*vblock*vblock)
cmp      call w_alloc(voa,vblock*vblock*nnoab(IUHF),'voaT3loopa')
      call GetMem('loopa_voa','Allo','Real',
     & voa,vblock*vblock*nnoab(IUHF))
cmp      call w_alloc(vob,vblock*vblock*nnoab(IUHF),'vobT3loopa')
      call GetMem('loopa_vob','Allo','Real',
     & vob,vblock*vblock*nnoab(IUHF))
cmp      call w_alloc(voc,vblock*vblock*nnoab(IUHF),'vocT3loopa')
      call GetMem('loopa_voc','Allo','Real',
     & voc,vblock*vblock*nnoab(IUHF))
C  this is necessary
C  prefactors currently a formal allocation
cmp ?      call w_alloc(mi,1,'miT3loopa')
cmp ?      call w_alloc(mij,1,'T3loopa')
      call GetMem('loopa_mi','Allo','Real',mi,1)
      call GetMem('loopa_mij','Allo','Real',mij,1)

cmp!!!      endif     ! energ - initialization
      aset=(nga-1)*vblock
      adim=min(vblock,nuab(isp)-aset)
      bset=(ngb-1)*vblock
      bdim=min(vblock,nuab(isp)-bset)
      cset=(ngc-1)*vblock
      cdim=min(vblock,nuab(isp)-cset)
C
C case1 nga=ngb=ngc
C if memory is available loops over i,j,k, in subloops can be grouped !!!
      if(nga.eq.ngc) then
C
cmp      call t3_bt_aaa(nug,g(ka),g(ka),g(ka),g(la),g(mi),g(mij),
      call t3_bt_aaa(nug,Work(ka),Work(ka),Work(ka),Work(la),
     &Work(mi),Work(mij),
     $adim,N,noab(isp),nuab(isp),nnoab(iuhf),lu,iasblock,nga,oeh,
     $oep(aset+1),enx,Work(voa),
     $t1a(noab(isp)*aset+1),t1b(noab(isp)*aset+1),Work(t3a),Work(t3b),
     &ifvo)
cmp     $oep(aset+1),enx,g(voa),
cmp     $t1a(noab(isp)*aset+1),t1b(noab(isp)*aset+1),g(t3a),g(t3b),ifvo)
C
      elseif(nga.eq.ngb)then
cmp      call t3_bt_aac(nug,g(ka),g(kb),g(kc),g(la),g(lc),g(mi),g(mij),
      call t3_bt_aac(nug,Work(ka),Work(kb),Work(kc),Work(la),Work(lc),
     &Work(mi),Work(mij),
     $adim,cdim,N,noab(isp),nuab(isp),nnoab(iuhf),lu,iasblock,nga,ngc,
     $oeh,oep(aset+1),oep(cset+1),enx,Work(voa),Work(voc),
     $t1a(noab(isp)*aset+1),t1b(noab(isp)*aset+1),
     $t1a(noab(isp)*cset+1),t1b(noab(isp)*cset+1),
     $Work(t3a),Work(t3b),ifvo)
cmp     $oeh,oep(aset+1),oep(cset+1),enx,g(voa),g(voc),
cmp     $t1a(noab(isp)*aset+1),t1b(noab(isp)*aset+1),
cmp     $t1a(noab(isp)*cset+1),t1b(noab(isp)*cset+1),
cmp     $g(t3a),g(t3b),ifvo)
      elseif(ngb.eq.ngc)then
cmp      call t3_bt_acc(nug,g(ka),g(kb),g(kc),g(la),g(lc),g(mi),g(mij),
      call t3_bt_acc(nug,Work(ka),Work(kb),Work(kc),Work(la),Work(lc),
     &Work(mi),Work(mij),
     $adim,cdim,N,noab(isp),nuab(isp),nnoab(iuhf),lu,iasblock,nga,ngc,
     $oeh,oep(aset+1),oep(cset+1),enx,Work(voa),Work(voc),
     $t1a(noab(isp)*aset+1),t1b(noab(isp)*aset+1),
     $t1a(noab(isp)*cset+1),t1b(noab(isp)*cset+1),
     $Work(t3a),Work(t3b),ifvo)
cmp     $oeh,oep(aset+1),oep(cset+1),enx,g(voa),g(voc),
cmp     $t1a(noab(isp)*aset+1),t1b(noab(isp)*aset+1),
cmp     $t1a(noab(isp)*cset+1),t1b(noab(isp)*cset+1),
cmp     $g(t3a),g(t3b),ifvo)
      else
cmp      call t3_bt_abc(nug,g(ka),g(kb),g(kc),g(la),g(lb),g(lc),g(mi),
cmp     $g(mij),adim,bdim,cdim,N,noab(isp),nuab(isp),nnoab(iuhf),lu,
      call t3_bt_abc(nug,Work(ka),Work(kb),Work(kc),Work(la),Work(lb),
     &Work(lc),Work(mi),
     $Work(mij),adim,bdim,cdim,N,noab(isp),nuab(isp),nnoab(iuhf),lu,
     $iasblock,nga,ngb,ngc,
     $oeh,oep(aset+1),oep(bset+1),oep(cset+1),enx,Work(voa),Work(vob),
     &Work(voc),
     $t1a(noab(isp)*aset+1),t1b(noab(isp)*aset+1),t1a(noab(isp)*bset+1),
     $t1b(noab(isp)*bset+1),t1a(noab(isp)*cset+1),t1b(noab(isp)*cset+1),
     $Work(t3a),Work(t3b),ifvo)
cmp     $oeh,oep(aset+1),oep(bset+1),oep(cset+1),enx,g(voa),g(vob),g(voc),
cmp     $t1a(noab(isp)*aset+1),t1b(noab(isp)*aset+1),t1a(noab(isp)*bset+1),
cmp     $t1b(noab(isp)*bset+1),t1a(noab(isp)*cset+1),t1b(noab(isp)*cset+1),
cmp     $g(t3a),g(t3b),ifvo)
      endif   ! cases
      energ(isp)=energ(isp)+enx
cmp     write (*,'(A,i5,x,3(i3,1x),f21.19)') 'Tsk, nga, ngb, ngc, inc = ',
cmp     & jjj,nga,ngb,ngc,enx
cmp
c321     continue
cmp        write (6,*)
cmp        write (6,*) 'deallocating arrays in t3loopa'
cmp        write (6,*)
      call GetMem('loopa_mij','Free','Real',mij,1)
      call GetMem('loopa_mi','Free','Real',mi,1)
      call GetMem('loopa_voc','Free','Real',
     & voc,vblock*vblock*nnoab(IUHF))
      call GetMem('loopa_vob','Free','Real',
     & vob,vblock*vblock*nnoab(IUHF))
      call GetMem('loopa_voa','Free','Real',
     & voa,vblock*vblock*nnoab(IUHF))
      call GetMem('loopa_t3b','Free','Real',t3b,vblock*vblock*vblock)
      call GetMem('loopa_t3a','Free','Real',t3a,vblock*vblock*vblock)
      call GetMem('loopa_lc','Free','Real',lc,nnoab(IUHF)*vblock*n)
      call GetMem('loopa_lb','Free','Real',lb,nnoab(IUHF)*vblock*n)
      call GetMem('loopa_la','Free','Real',la,nnoab(IUHF)*vblock*n)
      if(nug.ne.1)then
      call GetMem('loopa_kc','Free','Real',kc,noab(isp)*vblock*vblock*n)
      call GetMem('loopa_kb','Free','Real',kb,noab(isp)*vblock*vblock*n)
      endif
      call GetMem('loopa_ka','Free','Real',ka,noab(isp)*vblock*vblock*n)
cmp
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_logical(lastcall)
        call Unused_integer(jjj)
      end if
      end
