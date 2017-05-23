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
cmp      SUBROUTINE t3loopa(oeh,oep,t1a,t1b,g,nga,ngb,ngc,vblock,energ,
      SUBROUTINE t3loopa(oeh,oep,t1a,t1b,nga,ngb,ngc,vblock,energ,
     $ isp,LU,ifvo,lastcall,scored,jjj,enx)
C  implemented integer offsets, PV, 16 may 2004.
      IMPLICIT NONE
#include "ndisk.fh"
#include "WrkSpc.fh"
cmp      real*8 g(*),energ(*),oeh(*),oep(*),enx,t1a(*),t1b(*)
      real*8 energ(*),oeh(*),oep(*),enx,t1a(*),t1b(*)
      integer nug
      integer isp,vblock,n,lu(*),nga,ngb,ngc,adim,bdim,cdim
      INTEGER iasblock(3),aset,bset,cset
      CHARACTER ich*1
      logical ifvo,lastcall,scored
      INTEGER IOPT,NOAB,NNOAB,NUAB,NNUAB,IUHF
cmp
        integer jjj
cmp
      COMMON/UHF/NOAB(2),NNOAB(3),NUAB(2),NNUAB(3),ICH(3)
      COMMON/IOIND/IOPT(96)
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
        call Unused_logical(jjj)
      end if
      end

      subroutine t3_bt_aaa(nug,ka,kb,kc,la,mi,mij,adim,N,noab,
     $nuab,nnoab,lu,iasblock,nga,oeh,oep,enx,voa,t1a,t1b,t3a,t3b,ifvo)
      implicit none
      real*8 one,zero
      parameter (one=1.d0,zero=0.d0)
      integer a,b,c,aa,ab,bc,ac
      real*8 XX,YY,enx,den,dena,denb,denc
      integer ndim,adim,noab,nuab,i,j,k,nga,iasblock(3),lu(2),N
      integer ias,nga_offset,nug_offset,jk,ij,ik,ii,jj,nug,nnoab
      logical ifvo
      real*8  ka(adim*(adim-1)/2,N,*),kb(*),kc(*)
      real*8  la(N*adim,nnoab),mi(*),mij(*)
!!      real*8  lb(N,adim,nnoab),lc(N,adim,nnoab)
      real*8  t3a(adim*(adim-1)/2,adim),t3b(adim*(adim-1)/2,adim)
!!      real*8  mi(adim*(adim-1)/2,adim,noab),mij(*)
      real*8  voa(adim*(adim-1)/2,nnoab)
      real*8  t1a(noab,*),t1b(noab,*),oeh(noab),oep(adim)
C
      if(adim.eq.1)return
cmp       write(6,*)'enter_aaa enx,nga,',nga,enx
      ndim=adim*(adim-1)/2
      call zeroma(t3b,1,ndim*adim)
      nug_offset=iasblock(1)*nug*(nug+1)/2
      ias=iasblock(2)*(nga-1)+1
      call multi_readir(la,nnoab*adim*N,lu(2),ias)
      ias=iasblock(2)*nug+iasblock(3)*(nga*(nga+1)/2-1)+1
      call multi_readir(voa,nnoab*ndim,lu(2),ias)
C
      nga_offset=iasblock(1)*(nga*(nga+1)/2-1)+1
         do i=1,noab
         ias=(i-1)*nug_offset+nga_offset
         call multi_readir(ka(1,1,i),N*ndim,lu(1),ias)
         enddo
         do i=3,noab
         ii=(i-1)*N*ndim+1
         jk=0
         do j=2,i-1
         jj=(j-1)*N*ndim+1
         ij=(i-1)*(i-2)/2+j
         ik=(i-1)*(i-2)/2
         do k=1,j-1
         jk=jk+1
         ik=ik+1
C  K_ab^ir x L_rc^jk
         call DGEMM_('N','N',ndim,adim,N,one,ka(1,1,i),ndim,la(1,jk),N,
     $                     zero,t3a,ndim)
C  K_ab^kr x L_rc^ij
         call DGEMM_('N','N',ndim,adim,N,one,ka(1,1,k),ndim,la(1,ij),N,
     $                     one,t3a,ndim)
C  -K_ab^jr x L_rc^ik
         call DGEMM_('N','N',ndim,adim,N,-one,ka(1,1,j),ndim,la(1,ik),N,
     $                     one,t3a,ndim)
         den=oeh(i)+oeh(j)+oeh(k)
         do a=3,adim
         aa=(a-1)*(a-2)/2
         dena=den-oep(a)
         bc=0
         do b=2,a-1
         denb=dena-oep(b)
         ab=aa+b
         do c=1,b-1
         bc=bc+1
         denc=denb-oep(c)
         ac=aa+c
         xx=t3a(ab,c)+t3a(bc,a)-t3a(ac,b)
         yy=xx/denc
         enx=enx+yy*xx
         t3b(ab,c)=yy
         t3b(bc,a)=yy
         t3b(ac,b)=-yy
         enddo
         enddo
         enddo
C ccsd(T) part vvoo*t3
         call DGEMM_('N','N',1,adim,ndim,one,voa(1,ij),1,
     $  t3b,ndim,one,t1a(k,1),noab)
         call DGEMM_('N','N',1,adim,ndim,one,voa(1,jk),1,
     $  t3b,ndim,one,t1a(i,1),noab)
         call DGEMM_('N','N',1,adim,ndim,-one,voa(1,ik),1,
     $  t3b,ndim,one,t1a(j,1),noab)
C ccsd(T) part t2*t3
         if(ifvo) then
         call DGEMM_('N','N',1,adim,ndim,one,ka(1,i,j),1,
     $  t3b,ndim,one,t1b(k,1),noab)
         call DGEMM_('N','N',1,adim,ndim,one,ka(1,j,k),1,
     $  t3b,ndim,one,t1b(i,1),noab)
         call DGEMM_('N','N',1,adim,ndim,-one,ka(1,i,k),1,
     $  t3b,ndim,one,t1b(j,1),noab)
         endif
         enddo !k
         enddo !j
         enddo !i
         return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_real_array(kb)
        call Unused_real_array(kc)
        call Unused_real_array(mi)
        call Unused_real_array(mij)
        call Unused_integer(nuab)
      end if
         end

      subroutine t3_bt_aac(nug,ka,kb,kc,la,lc,mi,mij,adim,cdim,
     $N,noab,nuab,nnoab,lu,iasblock,nga,ngc,oeh,oepa,oepc,enx,voa,voc,
     $t1aa,t1ba,t1ac,t1bc,t3a,t3b,ifvo)
      implicit none
      real*8 one,zero,sumt3
      logical ifvo
      parameter (one=1.d0,zero=0.d0)
      integer nadim,adim,noab,nuab,i,j,k,nga,iasblock(3),lu(2),N
      integer ias,nga_offset,ngc_offset,nug_offset,jk,ij,ik,ncdim,cdim
      integer a,b,c,ab,nug,nnoab,ngc,iasci,iasai !,kaka
      real*8 dena,denb,denc,xx,yy
      real*8  ka(adim*(adim-1)/2,n,*),kc(adim*cdim,n,*),kb(*),
     $la(N*adim,nnoab),lc(N*cdim,nnoab)
!!      real*8  lb(N,adim,nnoab),lc(N,adim,nnoab)
      real*8  t3a(cdim,*),t3b(cdim,adim,*),t1ac(noab,*),t1bc(noab,*)
      real*8  mi(*),mij(*),voa(adim*(adim-1)/2,*),voc(adim*cdim,*)
      real*8  t1aa(noab,*),t1ba(noab,*),enx,oeh(noab),oepa(adim),den
      real*8  oepc(cdim)
C
      sumt3=0.d0
      if(adim.eq.1)return
      nadim=adim*(adim-1)/2
      ncdim=adim*cdim
      nug_offset=iasblock(1)*nug*(nug+1)/2
      ias=iasblock(2)*(nga-1)+1
      call multi_readir(la,nnoab*adim*N,lu(2),ias)
      ias=iasblock(2)*(ngc-1)+1
      call multi_readir(lc,nnoab*cdim*N,lu(2),ias)
C reads vvoo
      nga_offset=iasblock(3)*(nga*(nga-1)/2+nga-1)+1
      ngc_offset=iasblock(3)*(nga*(nga-1)/2+ngc-1)+1
      ias=iasblock(2)*nug+nga_offset
      call multi_readir(voa,nnoab*nadim,lu(2),ias)
      ias=iasblock(2)*nug+ngc_offset
      call multi_readir(voc,nnoab*ncdim,lu(2),ias)
C
      nga_offset=iasblock(1)*(nga*(nga-1)/2+nga-1)+1
      ngc_offset=iasblock(1)*(nga*(nga-1)/2+ngc-1)+1
      do i=1,noab
      iasci=(i-1)*nug_offset+ngc_offset
      call multi_readir(kc(1,1,i),N*ncdim,lu(1),iasci)
      iasai=(i-1)*nug_offset+nga_offset
      call multi_readir(ka(1,1,i),N*nadim,lu(1),iasai)
      enddo
         do i=3,noab
         jk=0
         do j=2,i-1
         ij=(i-1)*(i-2)/2+j
         ik=(i-1)*(i-2)/2
         do k=1,j-1
         jk=jk+1
         ik=ik+1
!!         write(6,'(9I5)')i,j,k,iasai,iasaj,iasak,iasci,iascj,iasck
C  K_ab^ir x L_rc^jk
      call DGEMM_('T','T',cdim,nadim,N,one, lc(1,jk),N,ka(1,1,i),nadim,
     $                     zero,t3a,cdim)
C  K_ab^kr x L_rc^ij
      call DGEMM_('T','T',cdim,nadim,N,one, lc(1,ij),N,ka(1,1,k),nadim,
     $                     one,t3a,cdim)
C  -K_ab^jr x L_rc^ik
      call DGEMM_('T','T',cdim,nadim,N,-one,lc(1,ik),N,ka(1,1,j),nadim,
     $                     one,t3a,cdim)
C  K_bc^ir x L_ra^jk
      call DGEMM_('N','N',ncdim,adim,N,one,kc(1,1,i),ncdim,la(1,jk),N,
     $                     zero,t3b,ncdim)
C  K_bc^kr x L_ra^ij
      call DGEMM_('N','N',ncdim,adim,N,one,kc(1,1,k),ncdim,la(1,ij),N,
     $                     one,t3b,ncdim)
C  -K_bc^jr x L_ra^ik
      call DGEMM_('N','N',ncdim,adim,N,-one,kc(1,1,j),ncdim,la(1,ik),N,
     $                     one,t3b,ncdim)
!!
       ab=0
      do a=2,adim
      do b=1,a-1
      ab=ab+1
      call daxpy_(cdim,-1.d0,t3b(1,a,b),1,t3a(1,ab),1)
      call daxpy_(cdim,1.d0,t3b(1,b,a),1,t3a(1,ab),1)
      enddo
      enddo
         den=oeh(i)+oeh(j)+oeh(k)
      ab=0
c!      kaka=0
      do a=2,adim
      dena=den-oepa(a)
      do b=1,a-1
      ab=ab+1
      denb=dena-oepa(b)
      do c=1,cdim
      denc=denb-oepc(c)
c!      ab=ab+1
c!      kaka=kaka+1
c!      xx=t3a(ab,1)
      xx=t3a(c,ab)
      yy=xx/denc
      sumt3=sumt3+xx
      enx=enx+yy*xx
c!      t3a(ab,1)=yy
      t3a(c,ab)=yy
c!      write (*,*) ab,t3a(ab,1)
c!      write (*,*) kaka,t3a(c,ab)
      enddo
      enddo
      enddo
      call expa2_uhf(t3a,cdim,adim,-1,t3b)
C ccsd(T) part vvoo*t3
         call DGEMM_('N','T',1,cdim,nadim, one,voa(1,ij),1,
     $  t3a,cdim,one,t1ac(k,1),noab)
         call DGEMM_('N','T',1,cdim,nadim, one,voa(1,jk),1,
     $  t3a,cdim,one,t1ac(i,1),noab)
         call DGEMM_('N','T',1,cdim,nadim,-one,voa(1,ik),1,
     $  t3a,cdim,one,t1ac(j,1),noab)
         call DGEMM_('N','N',1,adim,ncdim,-one,voc(1,ij),1,
     $  t3b,ncdim,one,t1aa(k,1),noab)
         call DGEMM_('N','N',1,adim,ncdim,-one,voc(1,jk),1,
     $  t3b,ncdim,one,t1aa(i,1),noab)
         call DGEMM_('N','N',1,adim,ncdim, one,voc(1,ik),1,
     $  t3b,ncdim,one,t1aa(j,1),noab)
C ccsd(T) part t2*t3
         if(ifvo) then
         call DGEMM_('N','T',1,cdim,nadim, one,ka(1,i,j),1,
     $  t3a,cdim,one,t1bc(k,1),noab)
         call DGEMM_('N','T',1,cdim,nadim, one,ka(1,j,k),1,
     $  t3a,cdim,one,t1bc(i,1),noab)
         call DGEMM_('N','T',1,cdim,nadim,-one,ka(1,i,k),1,
     $  t3a,cdim,one,t1bc(j,1),noab)
         call DGEMM_('N','N',1,adim,ncdim,-one,kc(1,i,j),1,
     $  t3b,ncdim,one,t1ba(k,1),noab)
         call DGEMM_('N','N',1,adim,ncdim,-one,kc(1,j,k),1,
     $  t3b,ncdim,one,t1ba(i,1),noab)
         call DGEMM_('N','N',1,adim,ncdim, one,kc(1,i,k),1,
     $  t3b,ncdim,one,t1ba(j,1),noab)
         endif
          enddo !k
         enddo !j
         enddo !i
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_real_array(kb)
        call Unused_real_array(mi)
        call Unused_real_array(mij)
        call Unused_integer(nuab)
      end if
      end

       subroutine t3_bt_abc(nug,ka,kb,kc,la,lb,lc,mi,mij,adim,bdim,
     $cdim,N,noab,nuab,nnoab,lu,iasblock,nga,ngb,ngc,oeh,oepa,oepb,
     $oepc,enx,voa,vob,voc,t1aa,t1ba,t1ab,t1bb,t1ac,t1bc,t3a,t3b,ifvo)
      implicit none
      real*8 one,zero,ddot_,enx,den,dena,denb,denc
      parameter (one=1.d0,zero=0.d0)
      logical ifvo
      integer nadim,adim,bdim,nbdim,cdim,ncdim,noab,nuab,i,j,k
      integer nga,ngb,ngc,N,iasblock(3),lu(2),jk,ij,ik,ac,ca,a,b,c
      integer ias,iasai,iasbi,iasci
c     integer iasaj,iasak,iascj,iasck
      integer nga_offset,ngb_offset,ngc_offset,nug_offset,nug,nnoab
      real*8  ka(adim*bdim,N,*),kb(bdim*cdim,N,*),kc(adim*cdim,N,*)
      real*8  la(N*adim,nnoab),lb(N*bdim,nnoab),lc(N*cdim,nnoab)
      real*8  t1aa(noab,*),t1ba(noab,*),t1ab(noab,*),t1bb(noab,*)
      real*8  mi(*),t1ac(noab,*),t1bc(noab,*)
      real*8  mij(*),voa(adim*bdim,*),vob(bdim*cdim,*),voc(adim*cdim,*)
      real*8  oeh(noab),oepa(adim),oepb(bdim),oepc(cdim),t3a(*),t3b(*)
C
!!       write(6,*)'enter_abc enx,nga,ngb,ngc',nga,ngb,ngc,enx
      nadim=adim*bdim
      nbdim=bdim*cdim
      ncdim=adim*cdim
      nug_offset=iasblock(1)*nug*(nug+1)/2
C reads la's
      ias=iasblock(2)*(nga-1)+1
      call multi_readir(la,nnoab*adim*N,lu(2),ias)
      ias=iasblock(2)*(ngb-1)+1
      call multi_readir(lb,nnoab*bdim*N,lu(2),ias)
      ias=iasblock(2)*(ngc-1)+1
      call multi_readir(lc,nnoab*cdim*N,lu(2),ias)
C reads vvoo
      nga_offset=iasblock(3)*(nga*(nga-1)/2+ngb-1)+1
      ngb_offset=iasblock(3)*(ngb*(ngb-1)/2+ngc-1)+1
      ngc_offset=iasblock(3)*(nga*(nga-1)/2+ngc-1)+1
      ias=iasblock(2)*nug+nga_offset
      call multi_readir(voa,nnoab*nadim,lu(2),ias)
      ias=iasblock(2)*nug+ngb_offset
      call multi_readir(vob,nnoab*nbdim,lu(2),ias)
      ias=iasblock(2)*nug+ngc_offset
      call multi_readir(voc,nnoab*ncdim,lu(2),ias)
C
      nga_offset=iasblock(1)*(nga*(nga-1)/2+ngb-1)+1
      ngb_offset=iasblock(1)*(ngb*(ngb-1)/2+ngc-1)+1
      ngc_offset=iasblock(1)*(nga*(nga-1)/2+ngc-1)+1
      do i=1,noab
      iasci=(i-1)*nug_offset+ngc_offset
      call multi_readir(kc(1,1,i),N*ncdim,lu(1),iasci)
      iasbi=(i-1)*nug_offset+ngb_offset
      call multi_readir(kb(1,1,i),N*nbdim,lu(1),iasbi)
      iasai=(i-1)*nug_offset+nga_offset
      call multi_readir(ka(1,1,i),N*nadim,lu(1),iasai)
      enddo
         do i=3,noab
         jk=0
         do j=2,i-1
         ij=(i-1)*(i-2)/2+j
         ik=(i-1)*(i-2)/2
         do k=1,j-1
         jk=jk+1
         ik=ik+1
C  K_ba^ir x L_rc^jk
       call DGEMM_('N','N',nadim,cdim,N,one,ka(1,1,i),nadim,lc(1,jk),N,
     $                     zero,t3a,nadim)
C  K_ba^kr x L_rc^ij
       call DGEMM_('N','N',nadim,cdim,N,one,ka(1,1,k),nadim,lc(1,ij),N,
     $                     one,t3a,nadim)
C  -K_ba^jr x L_rc^ik
       call DGEMM_('N','N',nadim,cdim,N,-one,ka(1,1,j),nadim,lc(1,ik),N,
     $                     one,t3a,nadim)
C  -K_ca^ir x L_rb^jk   ! in the matrix as b,c,a
       call DGEMM_('T','T',bdim,ncdim,N,-one,lb(1,jk),N,kc(1,1,i),ncdim,
     $                     zero,t3b,bdim)
C  -K_ca^kr x L_rb^ij
       call DGEMM_('T','T',bdim,ncdim,N,-one,lb(1,ij),N,kc(1,1,k),ncdim,
     $                     one ,t3b,bdim)
C  K_ca^jr x L_rb^ik
       call DGEMM_('T','T',bdim,ncdim,N,one,lb(1,ik),N,kc(1,1,j),ncdim,
     $                     one ,t3b,bdim)
!!
         ac=1-bdim
         do a=1,adim
         ca=(a-1)*bdim+1-nadim
         do c=1,cdim
         ac=ac+bdim
         ca=ca+nadim
         call daxpy_(bdim,1.d0,t3a(ca),1,t3b(ac),1)
         enddo
         enddo
C  K_cb^ir x L_ra^jk   ! in the matrix as b,c,a
       call DGEMM_('N','N',nbdim,adim,N,one,kb(1,1,i),nbdim,la(1,jk),N,
     $                                  zero,t3a,nbdim)
C  K_cb^kr x L_ra^ij
       call DGEMM_('N','N',nbdim,adim,N,one,kb(1,1,k),nbdim,la(1,ij),N,
     $                                  one ,t3a,nbdim)
C  -K_cb^jr x L_ra^ik
       call DGEMM_('N','N',nbdim,adim,N,-one,kb(1,1,j),nbdim,la(1,ik),N,
     $                                  one ,t3a,nbdim)
         den=oeh(i)+oeh(j)+oeh(k)
         ac=1
         do a=1,adim
         ca=(a-1)*nbdim+1
         do c=1,cdim
         call daxpy_(bdim,1.d0,t3a(ca),cdim,t3b(ac),1)
         ac=ac+bdim
         ca=ca+1
         enddo
         enddo
         call dcopy_(cdim*bdim*adim,t3b,1,t3a,1)
          ac=0
         do a=1,adim
         dena=den-oepa(a)
         do c=1,cdim
         denc=dena-oepc(c)
         do b=1,bdim
         denb=denc-oepb(b)
         ac=ac+1
         t3b(ac)=t3b(ac)/denb
         enddo
         enddo
         enddo
         enx=enx+ddot_(cdim*bdim*adim,t3b,1,t3a,1)
         call ex23(t3b,t3a,bdim,cdim,adim,1)
C t3a  bac
C blok1   voa ka
         call DGEMM_('N','N',1,cdim,nadim, one,voa(1,ij),1,
     $  t3a,nadim,one,t1ac(k,1),noab)
         call DGEMM_('N','N',1,cdim,nadim, one,voa(1,jk),1,
     $  t3a,nadim,one,t1ac(i,1),noab)
         call DGEMM_('N','N',1,cdim,nadim,-one,voa(1,ik),1,
     $  t3a,nadim,one,t1ac(j,1),noab)
C blok 2 t3b bca
            call DGEMM_('N','T',1,bdim,ncdim,-one, voc(1,ij),1,
     $  t3b,bdim,one,t1ab(k,1),noab)
            call DGEMM_('N','T',1,bdim,ncdim,-one, voc(1,jk),1,
     $  t3b,bdim,one,t1ab(i,1),noab)
            call DGEMM_('N','T',1,bdim,ncdim,one,voc(1,ik),1,
     $  t3b,bdim,one,t1ab(j,1),noab)
            if(ifvo)then
             call DGEMM_('N','N',1,cdim,nadim, one,ka(1,i,j),1,
     $  t3a,nadim,one,t1bc(k,1),noab)
             call DGEMM_('N','N',1,cdim,nadim, one,ka(1,j,k),1,
     $  t3a,nadim,one,t1bc(i,1),noab)
             call DGEMM_('N','N',1,cdim,nadim,-one,ka(1,i,k),1,
     $  t3a,nadim,one,t1bc(j,1),noab)
            call DGEMM_('N','T',1,bdim,ncdim,-one, kc(1,i,j),1,
     $  t3b,bdim,one,t1bb(k,1),noab)
            call DGEMM_('N','T',1,bdim,ncdim,-one, kc(1,j,k),1,
     $  t3b,bdim,one,t1bb(i,1),noab)
            call DGEMM_('N','T',1,bdim,ncdim, one,kc(1,i,k),1,
     $  t3b,bdim,one,t1bb(j,1),noab)
             endif
C part 3 acb in t3b
            call transm(t3a,t3b,bdim,ncdim)
              call DGEMM_('N','T',1,adim,nbdim,one, vob(1,ij),1,
     $  t3b,adim,one,t1aa(k,1),noab)
               call DGEMM_('N','T',1,adim,nbdim,one, vob(1,jk),1,
     $  t3b,adim,one,t1aa(i,1),noab)
               call DGEMM_('N','T',1,adim,nbdim,-one,vob(1,ik),1,
     $  t3b,adim,one,t1aa(j,1),noab)
               if(ifvo)then
               call DGEMM_('N','T',1,adim,nbdim,one,  kb(1,i,j),1,
     $  t3b,adim,one,t1ba(k,1),noab)
               call DGEMM_('N','T',1,adim,nbdim,one,  kb(1,j,k),1,
     $  t3b,adim,one,t1ba(i,1),noab)
               call DGEMM_('N','T',1,adim,nbdim,-one, kb(1,i,k),1,
     $  t3b,adim,one,t1ba(j,1),noab)
               endif
         enddo !k
         enddo !j
         enddo !i
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_real_array(mi)
        call Unused_real_array(mij)
        call Unused_integer(nuab)
      end if
      end

      subroutine t3_bt_acc(nug,ka,kb,kc,la,lc,mi,mij,adim,cdim,
     $N,noab,nuab,nnoab,lu,iasblock,nga,ngc,oeh,oepa,oepc,enx,voa,voc,
     $t1aa,t1ba,t1ac,t1bc,t3a,t3b,ifvo)
      implicit none
      real*8 one,zero,sumt3
      parameter (one=1.d0,zero=0.d0)
      logical ifvo
      integer nadim,adim,noab,nuab,i,j,k,nga,iasblock(3),lu(2),N
      integer ias,nga_offset,ngc_offset,nug_offset,jk,ij,ik,ncdim,cdim
      integer a,b,c,bc,abc,nug,nnoab,ngc,iasci,iasai
      real*8  ka(adim*cdim,N,*),kc(cdim*(cdim-1)/2,N,*),kb(*)
      real*8  t3a(cdim,cdim,*),t3b(*),t1ac(noab,*),t1bc(noab,*)
      real*8  lc(N*cdim,nnoab),la(N*adim,nnoab)
      real*8  mi(*),mij(*),voa(adim*cdim,*),voc(cdim*(cdim-1)/2,*)
      real*8  t1aa(noab,*),t1ba(noab,*),enx,oeh(noab),oepa(adim)
      real*8  oepc(cdim),den,dena,denb,denc,xx,yy
C
      sumt3=0.d0
      if(cdim.eq.1)return
      ncdim=cdim*(cdim-1)/2
      nadim=adim*cdim
      nug_offset=iasblock(1)*nug*(nug+1)/2
      ias=iasblock(2)*(nga-1)+1
      call multi_readir(la,nnoab*adim*N,lu(2),ias)
      ias=iasblock(2)*(ngc-1)+1
      call multi_readir(lc,nnoab*cdim*N,lu(2),ias)
C reads vvoo
      nga_offset=iasblock(3)*(nga*(nga-1)/2+ngc-1)+1
      ngc_offset=iasblock(3)*(ngc*(ngc-1)/2+ngc-1)+1
      ias=iasblock(2)*nug+nga_offset
      call multi_readir(voa,nnoab*nadim,lu(2),ias)
      ias=iasblock(2)*nug+ngc_offset
      call multi_readir(voc,nnoab*ncdim,lu(2),ias)
C
      nga_offset=iasblock(1)*(nga*(nga-1)/2+ngc-1)+1
      ngc_offset=iasblock(1)*(ngc*(ngc-1)/2+ngc-1)+1
      do i=1,noab
      iasci=(i-1)*nug_offset+ngc_offset
      call multi_readir(kc(1,1,i),N*ncdim,lu(1),iasci)
      iasai=(i-1)*nug_offset+nga_offset
      call multi_readir(ka(1,1,i),N*nadim,lu(1),iasai)
      enddo
      do i=3,noab
         jk=0
         do j=2,i-1
         ij=(i-1)*(i-2)/2+j
         ik=(i-1)*(i-2)/2
         do k=1,j-1
         jk=jk+1
         ik=ik+1
C  K_ab^ir x L_rc^jk   (stored as   c,b,a)
        call DGEMM_('T','T',cdim,nadim,N,one, lc(1,jk),N,
     $  ka(1,1,i),nadim,zero,t3a,cdim)
C  K_ab^kr x L_rc^ij
        call DGEMM_('T','T',cdim,nadim,N,one, lc(1,ij),N,ka(1,1,k),
     $              nadim,one,t3a,cdim)
C  -K_ab^jr x L_rc^ik
        call DGEMM_('T','T',cdim,nadim,N,-one,lc(1,ik),N,ka(1,1,j),
     $              nadim,one,t3a,cdim)
C
C  K_bc^ir x L_ra^jk
        call DGEMM_('N','N',ncdim,adim,N,one,kc(1,1,i),ncdim,la(1,jk),N,
     $                     zero,t3b,ncdim)
C  K_bc^kr x L_ra^ij
        call DGEMM_('N','N',ncdim,adim,N,one,kc(1,1,k),ncdim,la(1,ij),N,
     $                     one,t3b,ncdim)
C  -K_bc^jr x L_ra^ik
        call DGEMM_('N','N',ncdim,adim,N,-one,kc(1,1,j),ncdim,
     &        la(1,ik),N,
     $                     one,t3b,ncdim)
!!
         den=oeh(i)+oeh(j)+oeh(k)
         abc=0
         do a=1,adim
      dena=den-oepa(a)
       bc=0
      do b=2,cdim
      denb=dena-oepc(b)
      do c=1,b-1
      denc=denb-oepc(c)
      bc=bc+1
      abc=abc+1
      xx=t3b(abc)-t3a(b,c,a)+t3a(c,b,a)
      yy=xx/denc
      t3b(abc)=yy
      sumt3=sumt3+xx
      enx=enx+yy*xx
      enddo
      enddo
      enddo
      call expa1_uhf(t3b,adim,cdim,-1,t3a)
C t3a  bac
C blok1   voa ka
            call DGEMM_('N','T',1,cdim,nadim,-one,voa(1,ij),1,
     $  t3a,cdim,one,t1ac(k,1),noab)
            call DGEMM_('N','T',1,cdim,nadim,-one,voa(1,jk),1,
     $  t3a,cdim,one,t1ac(i,1),noab)
            call DGEMM_('N','T',1,cdim,nadim, one,voa(1,ik),1,
     $  t3a,cdim,one,t1ac(j,1),noab)
C blok 2 t3b bca
            call DGEMM_('N','N',1,adim,ncdim,one, voc(1,ij),1,
     $  t3b,ncdim,one,t1aa(k,1),noab)
            call DGEMM_('N','N',1,adim,ncdim,one, voc(1,jk),1,
     $  t3b,ncdim,one,t1aa(i,1),noab)
            call DGEMM_('N','N',1,adim,ncdim,-one,voc(1,ik),1,
     $  t3b,ncdim,one,t1aa(j,1),noab)
            if(ifvo)then
            call DGEMM_('N','T',1,cdim,nadim,-one,ka(1,i,j),1,
     $  t3a,cdim,one,t1bc(k,1),noab)
            call DGEMM_('N','T',1,cdim,nadim,-one,ka(1,j,k),1,
     $  t3a,cdim,one,t1bc(i,1),noab)
            call DGEMM_('N','T',1,cdim,nadim, one,ka(1,i,k),1,
     $  t3a,cdim,one,t1bc(j,1),noab)
            call DGEMM_('N','N',1,adim,ncdim,one, kc(1,i,j),1,
     $  t3b,ncdim,one,t1ba(k,1),noab)
            call DGEMM_('N','N',1,adim,ncdim,one, kc(1,j,k),1,
     $  t3b,ncdim,one,t1ba(i,1),noab)
            call DGEMM_('N','N',1,adim,ncdim,-one,kc(1,i,k),1,
     $  t3b,ncdim,one,t1ba(j,1),noab)
            endif
         enddo !k
         enddo !j
         enddo !i
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_real_array(kb)
        call Unused_real_array(mi)
        call Unused_real_array(mij)
        call Unused_integer(nuab)
      end if
      end
