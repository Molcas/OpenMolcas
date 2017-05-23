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
cmp      SUBROUTINE t3loopb(oeh,oep,t1a,t1b,g,nga,ngb,ngc,vblock,energ,
      SUBROUTINE t3loopb(oeh,oep,t1a,t1b,nga,ngb,ngc,vblock,energ,
     $ isp,LU,ifvo,lastcall,scored,jjj,enx)
C  implemented integer offsets, PV, 16 may 2004.
      IMPLICIT NONE
#include "ndisk.fh"
#include "WrkSpc.fh"
cmp      real*8 g(*),energ(*),oeh(*),oep(*),enx,t1a(*),t1b(*)
      real*8 energ(*),oeh(*),oep(*),enx,t1a(*),t1b(*)
      logical ifvo,lastcall,scored
      integer isp,vblock,n,lu(6),nga,ngb,ngc,adim,bdim,cdim
      integer en_offset_ah,en_offset_bh,en_offset_ap,en_offset_bp
      integer t1_offset_a,t1_offset_b
      integer nuga,nugc
      INTEGER iasblock(5),aset,bset,cset
      CHARACTER ich*1
      INTEGER IOPT,NOAB,NNOAB,NUAB,NNUAB
cmp
        integer jjj
cmp
      COMMON/UHF/NOAB(2),NNOAB(3),NUAB(2),NNUAB(3),ICH(3)
      COMMON/IOIND/IOPT(96)
      integer kab,kca,kcb,kac,kbc,kc,la,lb,lxa,lxb,lxc,t3a,t3b,
     $     vab,vbc,vac,mij,mi
cmp      SAVE      kab,kca,kcb,kac,kbc,kc,la,lb,lxa,lxb,lxc,t3a,t3b,
cmp     $     vab,vbc,vac,mij,mi,
cmp     $     iasblock,nuga,nugc
C
      en_offset_ah=(isp-1)*noab(1)
      en_offset_bh=(2-isp)*noab(1)
      en_offset_ap=(isp-1)*nuab(1)
      en_offset_bp=(2-isp)*nuab(1)
      t1_offset_a=(isp-1)*noab(1)*nuab(1)
      t1_offset_b=(2-isp)*noab(1)*nuab(1)
      N=noab(isp)+nuab(isp)
c
cmp
        scored=.true.
cmp!!!        if (lastcall) goto 321
cmp
      enx=0.d0
cmpn        write (*,*) 'Bef isp, energ(isp), enx = ',isp,energ(isp),enx
cmp!!!      if(energ(isp).eq.0.d0)then
C this is a first entry - initialization (makes no harm if reapeated)
      nuga=nuab(isp)/vblock
      if((nuga*vblock).lt.nuab(isp))nuga=nuga+1
!!      write(6,*)'first,nuga,vblock',nuga,vblock
      nugc=nuab(3-isp)/vblock
      if((nugc*vblock).lt.nuab(3-isp))nugc=nugc+1
!!      write(6,*)'first,nugc,vblock',nugc,vblock

      iasblock(1)=vblock*vblock*N/nblock
      if((iasblock(1)*nblock).lt.(vblock*vblock*N))
     $iasblock(1)=iasblock(1)+1

      iasblock(2)=nnoab(isp)*vblock*N/nblock
      if((iasblock(2)*nblock).lt.(nnoab(isp)*vblock*N))
     $iasblock(2)=iasblock(2)+1

      iasblock(3)=nnoab(3)*vblock*N/nblock
      if((iasblock(3)*nblock).lt.(nnoab(3)*vblock*N))
     $iasblock(3)=iasblock(3)+1

      iasblock(4)=nnoab(isp)*vblock*vblock/nblock
      if((iasblock(4)*nblock).lt.(nnoab(isp)*vblock*vblock))
     $iasblock(4)=iasblock(4)+1

      iasblock(5)=nnoab(3)*vblock*vblock/nblock
      if((iasblock(5)*nblock).lt.(nnoab(3)*vblock*vblock))
     $iasblock(5)=iasblock(5)+1

cmp      call w_rescope(G,'GT3loopb')
cmp      call w_free(g,0,'GT3loopb')
c  allocations
      if(nuga.ne.1)then
cmp      call w_alloc(kab,noab(isp)*vblock*vblock*n,'kaT3loopb')
        Call GetMem('loopb_kab','Allo','Real',kab,
     & noab(isp)*vblock*vblock*n)
cmp      call w_alloc(kcb,noab(isp)*vblock*vblock*n,'kbT3loopb')
        Call GetMem('loopb_kcb','Allo','Real',kcb,
     & noab(isp)*vblock*vblock*n)
cmp      call w_alloc(kbc,vblock*vblock*n,'kcT3loopb')
        Call GetMem('loopb_kbc','Allo','Real',kbc,
     & vblock*vblock*n)
      else
cmp      call w_alloc(kab,noab(isp)*N*nnuab(isp),'kaT3loopb')
        Call GetMem('loopb_kab','Allo','Real',kab,
     & noab(isp)*N*nnuab(isp))
      endif
cmp      call w_alloc(kac,vblock*vblock*n,'kcT3loopb')
        Call GetMem('loopb_kac','Allo','Real',kac,
     & vblock*vblock*n)
cmp      call w_alloc(kca,noab(isp)*vblock*vblock*n,'kbT3loopb')
        Call GetMem('loopb_kca','Allo','Real',kca,
     & noab(isp)*vblock*vblock*n)
cmp      call w_alloc(kc,vblock*vblock*n,'kcT3loopb')
        Call GetMem('loopb_kc','Allo','Real',kc,
     & vblock*vblock*n)
cmp      call w_alloc(la,nnoab(isp)*vblock*n,'laT3loopb')
        Call GetMem('loopb_la','Allo','Real',la,
     & nnoab(isp)*vblock*n)
cmp      call w_alloc(lxa,nnoab(3)*vblock*n,'lbaT3loopb')
        Call GetMem('loopb_lxa','Allo','Real',lxa,
     & nnoab(3)*vblock*n)
cmpn        write (*,*) 'check lxa'
cmpn        write (*,*) 'nnoab(1), (2), (3) = ',
cmpn     &               nnoab(1),nnoab(2),nnoab(3)
cmpn        write (*,*) 'nnuab(1), (2), (3) = ',
cmpn     &               nnuab(1),nnuab(2),nnuab(3)
       if(nuga.ne.1)then
cmp      call w_alloc(lb,nnoab(isp)*vblock*n,'lbT3loopb')
        Call GetMem('loopb_lb','Allo','Real',lb,
     & nnoab(isp)*vblock*n)
cmp      call w_alloc(lxb,nnoab(3)*vblock*n,'labT3loopb')
        Call GetMem('loopb_lxb','Allo','Real',lxb,
     & nnoab(3)*vblock*n)
       endif
cmp      call w_alloc(lxc,nnoab(3)*vblock*n,'lacT3loopb')
        Call GetMem('loopb_lxc','Allo','Real',lxc,
     & nnoab(3)*vblock*n)
cmp      call w_alloc(t3a,vblock*vblock*vblock,'t3aT3loopb')
        Call GetMem('loopb_t3a','Allo','Real',t3a,
     & vblock*vblock*vblock)
cmp      call w_alloc(t3b,vblock*vblock*vblock,'t3bT3loopb')
        Call GetMem('loopb_t3b','Allo','Real',t3b,
     & vblock*vblock*vblock)
cmp      call w_alloc(vac,vblock*vblock*nnoab(3),'vbcT3loopb')
        Call GetMem('loopb_vac','Allo','Real',vac,
     & vblock*vblock*nnoab(3))
      if(nuga.ne.1) then
cmp      call w_alloc(vab,vblock*vblock*nnoab(isp),'vabT3loopb')
        Call GetMem('loopb_vab','Allo','Real',vab,
     & vblock*vblock*nnoab(isp))
cmp      call w_alloc(vbc,vblock*vblock*nnoab(3),'vacT3loopb')
        Call GetMem('loopb_vbc','Allo','Real',vbc,
     & vblock*vblock*nnoab(3))
      else
cmp      call w_alloc(vab,nnoab(isp)*nnuab(isp),'vabT3loopb')
        Call GetMem('loopb_vab','Allo','Real',vab,
     & nnoab(isp)*nnuab(isp))
      endif
cmp      call w_alloc(mi,noab(isp)*vblock**3,'miT3loopb')
        Call GetMem('loopb_mi','Allo','Real',mi,
     & noab(isp)*vblock**3)
cmp
cmpn        write (*,*) 'mi check = ',noab(isp)*vblock**3
cmpn        write (*,*) 'nnoab(1), (2), (3) = ',
cmpn     &               nnoab(1),nnoab(2),nnoab(3)
cmpn        write (*,*) 'nnuab(1), (2), (3) = ',
cmpn     &               nnuab(1),nnuab(2),nnuab(3)
cmp
cmp      call w_alloc(mij,N*vblock,'mijT3loopb')
        Call GetMem('loopb_mij','Allo','Real',mij,
     & N*vblock)
cmp!!!      endif     ! energ = 0 - initialization
C
      aset=(nga-1)*vblock
      adim=min(vblock,nuab(isp)-aset)
      bset=(ngb-1)*vblock
      bdim=min(vblock,nuab(isp)-bset)
      cset=(ngc-1)*vblock
      cdim=min(vblock,nuab(3-isp)-cset)
C
C case1 nga=ngb=ngc
C
      if(nga.eq.ngb) then
C
C
cmp      call t3_bta_aac(nuga,nugc,g(kab),g(kca),g(kac),g(kc),g(la),g(lxa),
cmp     $g(lxc),g(mi),g(mij),adim,cdim,N,noab(isp),nuab(isp),
cmpn        write (*,*) 'ide call t3_bta_aac'
      call t3_bta_aac(nuga,nugc,Work(kab),Work(kca),Work(kac),Work(kc),
     &Work(la),Work(lxa),
     $Work(lxc),Work(mi),Work(mij),adim,cdim,N,noab(isp),nuab(isp),
     $noab(3-isp),nuab(3-isp),lu,iasblock,nga,ngc,
     $oeh(en_offset_ah+1),oeh(en_offset_bh+1),
     $oep(aset+en_offset_ap+1),oep(cset+en_offset_bp+1)
     $,enx,Work(vab),Work(vac),t1a(noab(isp)*aset+t1_offset_a+1),
     $t1b(noab(isp)*aset+t1_offset_a+1),
     $t1a(noab(3-isp)*cset+t1_offset_b+1),
     $t1b(noab(3-isp)*cset+t1_offset_b+1),Work(t3a),Work(t3b),ifvo)
cmp     $,enx,g(vab),g(vac),t1a(noab(isp)*aset+t1_offset_a+1),
cmp     $t1b(noab(isp)*aset+t1_offset_a+1),
cmp     $t1a(noab(3-isp)*cset+t1_offset_b+1),
cmp     $t1b(noab(3-isp)*cset+t1_offset_b+1),g(t3a),g(t3b),ifvo)
      else
cmp      call t3_bta_abc(nuga,nugc,g(kab),g(kcb),g(kca),g(kac),g(kbc),g(kc)
cmp     $,g(la),g(lb),
cmp     $g(lxa),g(lxb),g(lxc),g(mi),g(mij),adim,bdim,cdim,N,noab(isp),
cmpn        write (*,*) 'ide call t3_bta_abc'
cmp
      call t3_bta_abc(nuga,nugc,Work(kab),Work(kcb),Work(kca),Work(kac),
     &Work(kbc),Work(kc),Work(la),Work(lb),
     $Work(lxa),Work(lxb),Work(lxc),Work(mi),Work(mij),adim,bdim,cdim,
     &N,noab(isp),
     $nuab(isp),noab(3-isp),nuab(3-isp),lu,iasblock,nga,ngb,ngc,
     $oeh(en_offset_ah+1),oeh(en_offset_bh+1),oep(aset+en_offset_ap+1),
     $oep(bset+en_offset_ap+1),oep(cset+en_offset_bp+1)
     $,enx,Work(vab),Work(vbc),Work(vac),t1a(noab(isp)*aset+
     &t1_offset_a+1),
     $t1b(noab(isp)*aset+t1_offset_a+1),
     $t1a(noab(isp)*bset+t1_offset_a+1),
     $t1b(noab(isp)*bset+t1_offset_a+1),
     $t1a(noab(3-isp)*cset+t1_offset_b+1),
     $t1b(noab(3-isp)*cset+t1_offset_b+1),Work(t3a),Work(t3b),ifvo)
cmp     $,enx,g(vab),g(vbc),g(vac),t1a(noab(isp)*aset+t1_offset_a+1),
cmp     $t1b(noab(isp)*aset+t1_offset_a+1),
cmp     $t1a(noab(isp)*bset+t1_offset_a+1),
cmp     $t1b(noab(isp)*bset+t1_offset_a+1),
cmp     $t1a(noab(3-isp)*cset+t1_offset_b+1),
cmp     $t1b(noab(3-isp)*cset+t1_offset_b+1),g(t3a),g(t3b),ifvo)
      endif   ! cases
cmpn        write (*,*) 'isp, energ(isp), enx = ',isp,energ(isp),enx
      energ(isp)=energ(isp)+enx
cmp!!!        write (*,'(A,i5,x,3(i5,2x),f21.19)') 'Tsk, nga, ngb, ngc, inc = ',
cmp!!!     & jjj,nga,ngb,ngc,enx
c321     continue
cmp        write (6,*)
cmp        write (6,*) 'deallocating arrays in t3loob'
cmp        write (6,*)
cmp
        Call GetMem('loopb_mij','Free','Real',mij,
     & N*vblock)
        Call GetMem('loopb_mi','Free','Real',mi,
     & noab(isp)*vblock**3)
      if(nuga.ne.1) then
        Call GetMem('loopb_vbc','Free','Real',vbc,
     & vblock*vblock*nnoab(3))
        Call GetMem('loopb_vab','Free','Real',vab,
     & vblock*vblock*nnoab(isp))
      else
        Call GetMem('loopb_vab','Free','Real',vab,
     & nnoab(isp)*nnuab(isp))
      endif
        Call GetMem('loopb_vac','Free','Real',vac,
     & vblock*vblock*nnoab(3))
        Call GetMem('loopb_t3b','Free','Real',t3b,
     & vblock*vblock*vblock)
        Call GetMem('loopb_t3a','Free','Real',t3a,
     & vblock*vblock*vblock)
        Call GetMem('loopb_lxc','Free','Real',lxc,
     & nnoab(3)*vblock*n)
       if(nuga.ne.1)then
        Call GetMem('loopb_lxb','Free','Real',lxb,
     & nnoab(3)*vblock*n)
        Call GetMem('loopb_lb','Free','Real',lb,
     & nnoab(isp)*vblock*n)
       endif
        Call GetMem('loopb_lxa','Free','Real',lxa,
     & nnoab(3)*vblock*n)
        Call GetMem('loopb_la','Free','Real',la,
     & nnoab(isp)*vblock*n)
        Call GetMem('loopb_kc','Free','Real',kc,
     & vblock*vblock*n)
        Call GetMem('loopb_kca','Free','Real',kca,
     & noab(isp)*vblock*vblock*n)
        Call GetMem('loopb_kac','Free','Real',kac,
     & vblock*vblock*n)
      if(nuga.ne.1)then
        Call GetMem('loopb_kbc','Free','Real',kbc,
     & vblock*vblock*n)
        Call GetMem('loopb_kcb','Free','Real',kcb,
     & noab(isp)*vblock*vblock*n)
        Call GetMem('loopb_kab','Free','Real',kab,
     & noab(isp)*vblock*vblock*n)
      else
        Call GetMem('loopb_kab','Free','Real',kab,
     & noab(isp)*N*nnuab(isp))
      endif
cmp
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_logical(lastcall)
        call Unused_integer(jjj)
      end if
      end

      subroutine t3_bta_aac(nuga,nugc,kab,kca,kac,kc,la,lxa,lxc,mi,mij,
     $adim,cdim,N,noab_a,nuab_a,noab_b,nuab_b,lu,iasblock,nga,ngc,
     $oehi,oehk,oepa,oepc,enx,vab,vca,t1aa,t1ba,t1ac,t1bc,t3a,t3b,ifvo)
      implicit none
      real*8 one,zero,den,dena,denb,denc,enx,xx,yy
      parameter (one=1.d0,zero=0.d0)
      integer nadim,adim,ncdim,cdim,i,j,k,iasblock(5),lu(6),N
      integer noab_a,nuab_a,noab_b,nuab_b,nuga,nno_a,nnoab,nugc,
     $ngab_offset,ngca_offset,ngac_offset,nuga_offset,nugc_offset
      integer ias,jk,ij,ik,kj,ki,nga,ngc,a,b,c,ab,abb,bab
      integer iasabi,iascai,iasack
      real*8 kca(adim*cdim,N,*),kac(adim*cdim,N,*)
     $,kab(adim*(adim-1)/2,N,*)
     $,kc(*),la(N*adim,*),lxa(N*adim,*),lxc(N*cdim,*)
      real*8 t3a(*),t3b(*)
      real*8 mi(cdim*adim*(adim-1)/2,*),mij(*),
     $       vab(adim*(adim-1)/2,*),vca(adim*cdim,*)
      real*8 t1aa(noab_a,*),t1ba(noab_a,*),t1ac(noab_b,*),t1bc(noab_b,*)
      real*8 oehi(*),oehk(*),oepa(*),oepc(*)
      logical ifvo
C
C iasblock(1) > ka,kb,kc   iasblock(2) > la,lb iasblock(3) > lxa,lxc,lxb
      if(adim.eq.1)return
      nno_a=noab_a*(noab_a-1)/2
      nnoab=noab_a*noab_b
      nadim=adim*(adim-1)/2
      ncdim=adim*cdim
      nuga_offset=iasblock(1)*nuga*(nuga+1)/2
      nugc_offset=iasblock(1)*nuga*nugc
      ias=iasblock(2)*(nga-1)+1
      call multi_readir(la,nno_a*adim*N,lu(2),ias)
      ias=iasblock(3)*(nga-1)+1
      call multi_readir(lxa,nnoab*adim*N,lu(5),ias)
      ias=iasblock(3)*(ngc-1)+1
      call multi_readir(lxc,nnoab*cdim*N,lu(6),ias)
C vvoo ints reading
      ngab_offset=iasblock(4)*(nga*(nga-1)/2+nga-1)+1
      ias=iasblock(2)*nuga+ngab_offset
      call multi_readir(vab,nno_a*nadim,lu(2),ias)
      ngca_offset=iasblock(5)*(nugc*(nga-1)+ngc-1)+1
      ias=iasblock(2)*nuga+iasblock(4)*nuga*(nuga+1)/2+ngca_offset
      call multi_readir(vca,nnoab*adim*cdim,lu(2),ias)
!!      ngac_offset=iasblock(5)*(nuga*(ngc-1)+nga-1)+1
C end readin vvoo ints
      ngab_offset=iasblock(1)*(nga*(nga-1)/2+nga-1)+1
      ngac_offset=iasblock(1)*(nuga*(ngc-1)+nga-1)+1
      ngca_offset=iasblock(1)*(nugc*(nga-1)+ngc-1)+1

                do i=1,noab_a
                iasabi=(i-1)*nuga_offset+ngab_offset
                call multi_readir(kab(1,1,i),N*nadim,lu(1),iasabi)
                enddo
                do i=1,noab_a
                iascai=(i-1)*nugc_offset+ngca_offset
                call multi_readir(kca(1,1,i),N*ncdim,lu(3),iascai)
                enddo
         do k=1,noab_b
         do i=1,noab_a
         ik=(k-1)*noab_a +i
         ki=(i-1)*noab_b +k
C  K_ab^ir x L_rc^ik     cba
      call DGEMM_('T','T',cdim,nadim,N,one,lxc(1,ik),N,kab(1,1,i),nadim,
     $                      zero,mi(1,i),cdim)
C
C  K_ac^ir x L_rb^ki     cab
      call DGEMM_('N','N',ncdim,adim,N,one,kca(1,1,i),ncdim,lxa(1,ki),N,
     $                    zero,t3b,ncdim)
        ab=1
      do a=2,adim
      abb=(a-1)*cdim+1
      bab=(a-1)*ncdim+1
      do b=1,a-1
      call daxpy_(cdim,-1.d0,t3b(abb),1,mi(ab,i),1)
      call daxpy_(cdim,1.d0,t3b(bab),1,mi(ab,i),1)
      ab=ab+cdim
      abb=abb+ncdim
      bab=bab+cdim
      enddo
      enddo
          enddo    ! i
! end prefactors
         iasack=(k-1)*nugc_offset+ngac_offset
         call multi_readir(kac,N*ncdim,lu(4),iasack)
         ij=0
         do i=2,noab_a
         ki=(i-1)*noab_b +k
         ik=(k-1)*noab_a +i
         kj=k-noab_b
         jk=(k-1)*noab_a
         do j=1,i-1
         ij=ij+1
         kj=kj+noab_b
         jk=jk+1
C  K_bc^kr x L_ra^ij
         call DGEMM_('N','N',ncdim,adim,N,one,kac,ncdim,la(1,ij),N,
     $                     zero,t3a,ncdim)
C transpose the first two inicesd
         ab=1
         do a=1,adim
         call transm(t3a(ab),t3b(ab),adim,cdim)
         ab=ab+ncdim
         enddo
C  K_ab^ir x L_rc^jk -K_ab^jr x L_rc^ik
         call vsub(kab(1,1,j),1,kab(1,1,i),1,kc,1,N*nadim)
         call vadd(lxc(1,jk),1,lxc(1,ik),1,mij,1,N*cdim)
         call DGEMM_('T','T',cdim,nadim,N,one,mij,N,kc,nadim,
     $                     zero,t3a,cdim)
C  K_ab^ir x L_rc^jk
!!         call DGEMM_('T','T',cdim,nadim,N,one,lxc(1,jk),N,ka,nadim,
!!     $                     zero,t3a,cdim)
C  -K_ab^jr x L_rc^ik
!!         call DGEMM_('T','T',cdim,nadim,N,-one,lxc(1,ik),N,kb,nadim,
!!     $                      one,t3a,cdim)
C
C  K_bc^ir x L_ra^kj -K_bc^jr x L_ra^ki
         call vsub(kca(1,1,j),1,kca(1,1,i),1,kc,1,N*ncdim)
         call vadd(lxa(1,kj),1,lxa(1,ki),1,mij,1,N*adim)
         call DGEMM_('N','N',ncdim,adim,N,one,kc,ncdim,mij,N,
     $                     one,t3b,ncdim)
C  K_bc^ir x L_ra^kj
!!         call DGEMM_('N','N',ncdim,adim,N,one,ka,ncdim,lxa(1,kj),N,
!!     $                     one,t3b,ncdim)
C  -K_bc^jr x L_ra^ki
!!         call DGEMM_('N','N',ncdim,adim,N,-one,kb,ncdim,lxa(1,ki),N,
!!     $                     one,t3b,ncdim)
       ab=1
      do a=2,adim
      abb=(a-1)*cdim+1
      bab=(a-1)*ncdim+1
      do b=1,a-1
      call daxpy_(cdim,-1.d0,t3b(abb),1,t3a(ab),1)
      call daxpy_(cdim,1.d0,t3b(bab),1,t3a(ab),1)
      ab=ab+cdim
      abb=abb+ncdim
      bab=bab+cdim
      enddo
      enddo
!!
      call daxpy_(nadim*cdim,-1.d0,mi(1,i),1,t3a,1)
      call daxpy_(nadim*cdim,1.d0,mi(1,j),1,t3a,1)
         den=oehi(i)+oehi(j)+oehk(k)
      ab=0
      do a=2,adim
      dena=den-oepa(a)
      do b=1,a-1
      denb=dena-oepa(b)
      do c=1,cdim
      denc=denb-oepc(c)
      ab=ab+1
      xx=t3a(ab)
      yy=xx/denc
      enx=enx+yy*xx
      t3a(ab)=yy
!!      t1aa(j,a)=t1aa(j,a)-yy*vca((a-1)*cdim+c,ki)
      enddo
      enddo
      enddo
      call expa2_uhf(t3a,cdim,adim,-1,t3b)
         call DGEMM_('N','T',1,cdim,nadim, one,vab(1,ij),1,
     $  t3a,cdim,one,t1ac(k,1),noab_b)
         call DGEMM_('N','N',1,adim,ncdim, one,vca(1,kj),1,
     $  t3b,ncdim,one,t1aa(i,1),noab_a)
         call DGEMM_('N','N',1,adim,ncdim,-one,vca(1,ki),1,
     $  t3b,ncdim,one,t1aa(j,1),noab_a)
C ccsd(T) part t2*t3
         if(ifvo) then
         call DGEMM_('N','T',1,cdim,nadim, one,kab(1,i,j),1,
     $  t3a,cdim,one,t1bc(k,1),noab_b)
         call DGEMM_('N','N',1,adim,ncdim,-one,kca(1,k,j),1,
     $  t3b,ncdim,one,t1ba(i,1),noab_a)
         call DGEMM_('N','N',1,adim,ncdim, one,kca(1,k,i),1,
     $  t3b,ncdim,one,t1ba(j,1),noab_a)
         endif
         enddo !j
         enddo !i
         enddo !k
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(nuab_a)
        call Unused_integer(nuab_b)
      end if
      end

      subroutine t3_bta_abc(nuga,nugc,kab,kcb,kca,kac,kbc,kc,la,lb,lxa,
     $lxb,lxc,mi,mij,adim,bdim,cdim,N,noab_a,nuab_a,noab_b,nuab_b,lu,
     $iasblock,nga,ngb,ngc,oehi,oehk,oepa,oepb,oepc,enx,vab,vcb,vca,
     $t1aa,t1ba,t1ab,t1bb,t1ac,t1bc,t3a,t3b,ifvo)
      implicit none
cmp
c     integer imm
cmp
      real*8 one,zero,den,dena,denb,denc,enx,xx,yy
c     real*8 sumt3
      parameter (one=1.d0,zero=0.d0)
      integer nadim,adim,ncdim,cdim,bdim,nbdim,i,j,k,iasblock(5),lu(6),N
      integer noab_a,nuab_a,noab_b,nuab_b,nuga,nno_a,nnoab,nugc,
     $ngab_offset,ngca_offset,ngac_offset,nuga_offset,nugc_offset,
     $ngcb_offset,ngbc_offset
      integer ias,jk,ij,ik,kj,ki,nga,ngb,ngc,a,b,c,ab,ba
      integer iasabi,iascai,iasack,iascbi,iasbck
      real*8 kab(adim*bdim,N,*),kcb(cdim*bdim,N,*),kca(cdim*adim,N,*)
     $,kac(cdim*adim,N),kbc(cdim*bdim,N),kc(*),lxb(N*bdim,*)
      real*8 la(N*adim,*),lb(N*bdim,*),lxa(N*adim,*),lxc(N*cdim,*),
     $t3a(*),t3b(*),vab(adim*bdim,*),vca(adim*cdim,*),vcb(bdim*cdim,*)
      real*8 mi(adim*bdim*cdim,*),mij(*)
      real*8 t1aa(noab_a,*),t1ba(noab_a,*),t1ac(noab_b,*)
     $        ,t1bc(noab_b,*),t1ab(noab_a,*),t1bb(noab_a,*)
      real*8 oehi(*),oehk(*),oepa(*),oepb(*),oepc(*)
      logical ifvo
cmp
c
C iasblock(1) > ka,kb,kc   iasblock(2) > la,lb iasblock(3) > lxa,lxc,lxb
!!      sumt3=0.d0
      nno_a=noab_a*(noab_a-1)/2
      nnoab=noab_a*noab_b
      nadim=adim*bdim
      nbdim=bdim*cdim
      ncdim=adim*cdim
      nuga_offset=iasblock(1)*nuga*(nuga+1)/2
      nugc_offset=iasblock(1)*nuga*nugc
      ias=iasblock(2)*(nga-1)+1
      call multi_readir(la,nno_a*adim*N,lu(2),ias)
      ias=iasblock(2)*(ngb-1)+1
      call multi_readir(lb,nno_a*bdim*N,lu(2),ias)
      ias=iasblock(3)*(nga-1)+1
      call multi_readir(lxa,nnoab*adim*N,lu(5),ias)
      ias=iasblock(3)*(ngb-1)+1
      call multi_readir(lxb,nnoab*bdim*N,lu(5),ias)
      ias=iasblock(3)*(ngc-1)+1
      call multi_readir(lxc,nnoab*cdim*N,lu(6),ias)
C
C vvoo ints reading
      ngab_offset=iasblock(4)*(nga*(nga-1)/2+ngb-1)+1
      ias=iasblock(2)*nuga+ngab_offset
      call multi_readir(vab,nno_a*nadim,lu(2),ias)
      ngca_offset=iasblock(5)*(nugc*(nga-1)+ngc-1)+1
      ias=iasblock(2)*nuga+iasblock(4)*nuga*(nuga+1)/2+ngca_offset
      call multi_readir(vca,nnoab*ncdim,lu(2),ias)
      ngcb_offset=iasblock(5)*(nugc*(ngb-1)+ngc-1)+1
      ias=iasblock(2)*nuga+iasblock(4)*nuga*(nuga+1)/2+ngcb_offset
      call multi_readir(vcb,nnoab*nbdim,lu(2),ias)
C end readin vvoo ints
C
      ngab_offset=iasblock(1)*(nga*(nga-1)/2+ngb-1)+1
      ngac_offset=iasblock(1)*(nuga*(ngc-1)+nga-1)+1
      ngbc_offset=iasblock(1)*(nuga*(ngc-1)+ngb-1)+1
      ngca_offset=iasblock(1)*(nugc*(nga-1)+ngc-1)+1
      ngcb_offset=iasblock(1)*(nugc*(ngb-1)+ngc-1)+1
C saves reading:
      do i=1,noab_a
      iasabi=(i-1)*nuga_offset+ngab_offset
      call multi_readir(kab(1,1,i),N*nadim,lu(1),iasabi)
      enddo
      do i=1,noab_a
      iascai=(i-1)*nugc_offset+ngca_offset
      call multi_readir(kca(1,1,i),N*ncdim,lu(3),iascai)
cmp
cmp        write (*,*) 'ze tak teraz ju dam',i
cmp        call check_mat(kca(1,1,i),1,N*ncdim)
cmp
      enddo
      do i=1,noab_a
         iascbi=(i-1)*nugc_offset+ngcb_offset
         call multi_readir(kcb(1,1,i),N*nbdim,lu(3),iascbi)
      enddo
         do k=1,noab_b
         iasbck=(k-1)*nugc_offset+ngbc_offset
         call multi_readir(kbc,N*nbdim,lu(4),iasbck)
         iasack=(k-1)*nugc_offset+ngac_offset
         call multi_readir(kac,N*ncdim,lu(4),iasack)
C start calculating prefactors:
cmp
         do i=1,noab_a
         ik=(k-1)*noab_a +i
         ki=(i-1)*noab_b +k
cmp
C  K_ab^ir x L_rc^ik     cba
       call DGEMM_('T','T',cdim,nadim,N,one,lxc(1,ik),N,
     &  kab(1,1,i),nadim,
     $                      zero,mi(1,i),cdim)
C  K_bc^ir x L_ra^ki          cba
       call DGEMM_('N','N',nbdim,adim,N,one,kcb(1,1,i),nbdim,
     &  lxa(1,ki),N,
     $                     one ,mi(1,i),nbdim)
cmp
cmp        do imm=0,adim*nbdim-1
cmp        if (abs(mi(1+imm,i)).gt.10000) then
cmp          write (*,*) 'uz mi dojebane 2',imm+1,i,mi(1+imm,i)
cmp          stop
cmp        end if
cmp        end do
cmp
C
C  K_ac^ir x L_rb^ki     cab
      call DGEMM_('N','N',ncdim,bdim,N,one,kca(1,1,i),ncdim,
     &    lxb(1,ki),N,
     $                    zero,t3b,ncdim)
cmp
         ab=1
      do a=1,adim
      ba=(a-1)*cdim+1
      do b=1,bdim
      call daxpy_(cdim,-1.d0,t3b(ba),1,mi(ab,i),1)
      ab=ab+cdim
      ba=ba+ncdim
      enddo
      enddo
          enddo    ! i
cmp
! end prefactors
         ij=0
         do i=2,noab_a
         ki=(i-1)*noab_b +k
         ik=(k-1)*noab_a +i
         kj=k-noab_b
         jk=(k-1)*noab_a
         do j=1,i-1
         ij=ij+1
         kj=kj+noab_b
         jk=jk+1
cmp
C  K_ac^kr x L_rb^ij        bac
         call DGEMM_('T','T',bdim,ncdim,N,-one,lb(1,ij),N,kac,ncdim,
     $                     zero,t3b,bdim)
C  K_bc^kr x L_ra^ij        abc
         call DGEMM_('T','T',adim,nbdim,N,one,la(1,ij),N,kbc,nbdim,
     $                     zero,t3a,adim)
C transpose the first two indices
cmp
         ab=1
         do c=1,cdim
         ba=(c-1)*nadim+1
         do b=1,bdim
         call daxpy_(adim,1.d0,t3a(ab),1,t3b(ba),bdim)
         ba=ba+1
         ab=ab+adim
         enddo
         enddo
C t3b  bac
         call transm(t3b,t3a,nadim,cdim)
C      cba in t3a
C K_ab^ir x L_rc^jk  -K_ab^jr x L_rc^ik
         call vsub(kab(1,1,j),1,kab(1,1,i),1,kc,1,N*nadim)
C         call daxpy_(N*nadim,-one,kb,1,ka,1)
         call vadd(lxc(1,jk),1,lxc(1,ik),1,mij,1,N*cdim)
         call DGEMM_('T','T',cdim,nadim,N,one,mij,N,kc,nadim,
     $                      one,t3a,cdim)
C  K_ab^ir x L_rc^jk
!!         call DGEMM_('T','T',cdim,nadim,N,one,lxc(1,jk),N,ka,nadim,
!!     $                      one,t3a,cdim)
C  -K_ab^jr x L_rc^ik
!!         call DGEMM_('T','T',cdim,nadim,N,-one,lxc(1,ik),N,kb,nadim,
!!     $                      one,t3a,cdim)
C
C  K_bc^ir x L_ra^kj -K_bc^jr x L_ra^ki         cba
         call vsub(kcb(1,1,j),1,kcb(1,1,i),1,kc,1,N*nbdim)
         call vadd(lxa(1,kj),1,lxa(1,ki),1,mij,1,N*adim)
         call DGEMM_('N','N',nbdim,adim,N,one,kc,nbdim,mij,N,
     $                     one,t3a,nbdim)
C  K_bc^ir x L_ra^kj          cba
!!         call DGEMM_('N','N',nbdim,adim,N,one,ka,nbdim,lxa(1,kj),N,
!!     $                     one,t3a,nbdim)
C  -K_bc^jr x L_ra^ki         cba
!!         call DGEMM_('N','N',nbdim,adim,N,-one,kb,nbdim,lxa(1,ki),N,
!!     $                     one,t3a,nbdim)
C  K_ac^ir x L_rb^kj  -K_ac^jr x L_rb^ki      cab
         call vsub(kca(1,1,j),1,kca(1,1,i),1,kc,1,N*ncdim)
         call vadd(lxb(1,kj),1,lxb(1,ki),1,mij,1,N*bdim)
         call DGEMM_('N','N',ncdim,bdim,N,one,kc,ncdim,mij,N,
     $                    zero,t3b,ncdim)
C  K_ac^ir x L_rb^kj      cab
!!         call DGEMM_('N','N',ncdim,bdim,N,one,ka,ncdim,lxb(1,kj),N,
!!     $                    zero,t3b,ncdim)
C  -K_ac^jr x L_rb^ki     cab
!!         call DGEMM_('N','N',ncdim,bdim,N,-one,kb,ncdim,lxb(1,ki),N,
!!     $                     one,t3b,ncdim)
cmp
         ab=1
      do a=1,adim
      ba=(a-1)*cdim+1
      do b=1,bdim
      call daxpy_(cdim,-1.d0,t3b(ba),1,t3a(ab),1)
      ab=ab+cdim
      ba=ba+ncdim
      enddo
      enddo
c
cmp        do a=1,adim
cmp        do b=1,bdim
cmp        do c=1,cdim
cmp        ab=c+(b-1)*cdim+(a-1)*bdim*cdim
cmp        ba=c+(a-1)*cdim+(b-1)*adim*cdim
cmp        t3a(ab)=t3a(ab)-t3b(ba)
cmp        end do
cmp        end do
cmp        end do
c
cmp
c
      call daxpy_(nadim*cdim,-1.d0,mi(1,i),1,t3a,1)
      call daxpy_(nadim*cdim,1.d0,mi(1,j),1,t3a,1)
         den=oehi(i)+oehi(j)+oehk(k)
      ab=0
      do a=1,adim
      ba=ab+1
      dena=den-oepa(a)
      do b=1,bdim
      denb=dena-oepb(b)
      do c=1,cdim
      denc=denb-oepc(c)
      ab=ab+1
cmp        if ((i.eq.14).and.(j.eq.1).and.(k.eq.1)) then
cmp          write (*,'(A,3(i5,x),3(f18.10,x))')
cmp     &               'a,b,c, oepa(a), oepb(b), oepc(c)',
cmp     &                a,b,c,oepa(a),oepb(b),oepc(c)
cmp          write (*,*) 'denc = ',denc
cmp          write (*,*) 'ab, t3a(ab) ',ab,t3a(ab)
cmp        end if
      xx=t3a(ab)
      yy=xx/denc
!!      sumt3=sumt3+xx
      enx=enx+yy*xx
cmp        if ((i.eq.14).and.(j.eq.1).and.(k.eq.1)) then
cmp          write (*,'(A,3(f18.10,x))') 'xx,yy,enx = ',xx,yy,enx
cmp        end if
      t3a(ab)=yy
      enddo
      enddo
      call transm(t3a(ba),t3b(ba),cdim,bdim)
      enddo
      call DGEMM_('N','T',1,cdim,nadim, one,vab(1,ij),1,
     $  t3a,cdim,one,t1ac(k,1),noab_b)
       call DGEMM_('N','N',1,adim,nbdim,-one,vcb(1,kj),1,
     $  t3a,nbdim,one,t1aa(i,1),noab_a)
       call DGEMM_('N','N',1,adim,nbdim,one,vcb(1,ki),1,
     $  t3a,nbdim,one,t1aa(j,1),noab_a)
      call DGEMM_('N','T',1,bdim,ncdim,-one,vca(1,ki),1,
     $  t3b,bdim,one,t1ab(j,1),noab_a)
      call DGEMM_('N','T',1,bdim,ncdim,one,vca(1,kj),1,
     $  t3b,bdim,one,t1ab(i,1),noab_a)
       if(ifvo) then
      call DGEMM_('N','T',1,cdim,nadim, one,kab(1,i,j),1,
     $  t3a,cdim,one,t1bc(k,1),noab_b)
       call DGEMM_('N','N',1,adim,nbdim,one,kcb(1,k,j),1,
     $  t3a,nbdim,one,t1ba(i,1),noab_a)
       call DGEMM_('N','N',1,adim,nbdim,-one,kcb(1,k,i),1,
     $  t3a,nbdim,one,t1ba(j,1),noab_a)
      call DGEMM_('N','T',1,bdim,ncdim,one,kca(1,k,i),1,
     $  t3b,bdim,one,t1bb(j,1),noab_a)
      call DGEMM_('N','T',1,bdim,ncdim,-one,kca(1,k,j),1,
     $  t3b,bdim,one,t1bb(i,1),noab_a)
      endif
cmp        write (*,'(3(i5,x),f18.10)') i,j,k,enx
         enddo !j
         enddo !i
         enddo !k
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(nuab_a)
        call Unused_integer(nuab_b)
      end if
      end
c
c -----------
c
        subroutine check_mat(mat,dima,dimb)
c
        implicit none
        integer dima,dimb,i,j
        real*8 mat(dima,dimb)
c
        do i=1,dima
        do j=1,dimb
        if (abs(mat(i,j)).gt.10000) then
          write (6,*) 'i,j,mat(i,j) ',i,j,mat(i,j)
        end if
        end do
        end do
c
        return
        end
