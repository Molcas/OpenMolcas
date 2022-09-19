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
       subroutine t3_bt_abc(nug,ka,kb,kc,la,lb,lc,mi,mij,adim,bdim,     &
     &cdim,N,noab,nuab,nnoab,lu,iasblock,nga,ngb,ngc,oeh,oepa,oepb,     &
     &oepc,enx,voa,vob,voc,t1aa,t1ba,t1ab,t1bb,t1ac,t1bc,t3a,t3b,ifvo)
      implicit none
      real*8 one,zero,ddot_,enx,den,dena,denb,denc
      parameter (one=1.d0,zero=0.d0)
      logical ifvo
      integer nadim,adim,bdim,nbdim,cdim,ncdim,noab,nuab,i,j,k
      integer nga,ngb,ngc,N,iasblock(3),lu(2),jk,ij,ik,ac,ca,a,b,c
      integer ias,iasai,iasbi,iasci
!     integer iasaj,iasak,iascj,iasck
      integer nga_offset,ngb_offset,ngc_offset,nug_offset,nug,nnoab
      real*8  ka(adim*bdim,N,*),kb(bdim*cdim,N,*),kc(adim*cdim,N,*)
      real*8  la(N*adim,nnoab),lb(N*bdim,nnoab),lc(N*cdim,nnoab)
      real*8  t1aa(noab,*),t1ba(noab,*),t1ab(noab,*),t1bb(noab,*)
      real*8  mi(*),t1ac(noab,*),t1bc(noab,*)
      real*8  mij(*),voa(adim*bdim,*),vob(bdim*cdim,*),voc(adim*cdim,*)
      real*8  oeh(noab),oepa(adim),oepb(bdim),oepc(cdim),t3a(*),t3b(*)
!
!!       write(6,*)'enter_abc enx,nga,ngb,ngc',nga,ngb,ngc,enx
      nadim=adim*bdim
      nbdim=bdim*cdim
      ncdim=adim*cdim
      nug_offset=iasblock(1)*nug*(nug+1)/2
! reads la's
      ias=iasblock(2)*(nga-1)+1
      call multi_readir(la,nnoab*adim*N,lu(2),ias)
      ias=iasblock(2)*(ngb-1)+1
      call multi_readir(lb,nnoab*bdim*N,lu(2),ias)
      ias=iasblock(2)*(ngc-1)+1
      call multi_readir(lc,nnoab*cdim*N,lu(2),ias)
! reads vvoo
      nga_offset=iasblock(3)*(nga*(nga-1)/2+ngb-1)+1
      ngb_offset=iasblock(3)*(ngb*(ngb-1)/2+ngc-1)+1
      ngc_offset=iasblock(3)*(nga*(nga-1)/2+ngc-1)+1
      ias=iasblock(2)*nug+nga_offset
      call multi_readir(voa,nnoab*nadim,lu(2),ias)
      ias=iasblock(2)*nug+ngb_offset
      call multi_readir(vob,nnoab*nbdim,lu(2),ias)
      ias=iasblock(2)*nug+ngc_offset
      call multi_readir(voc,nnoab*ncdim,lu(2),ias)
!
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
!  K_ba^ir x L_rc^jk
       call DGEMM_('N','N',nadim,cdim,N,one,ka(1,1,i),nadim,lc(1,jk),N, &
     &                     zero,t3a,nadim)
!  K_ba^kr x L_rc^ij
       call DGEMM_('N','N',nadim,cdim,N,one,ka(1,1,k),nadim,lc(1,ij),N, &
     &                     one,t3a,nadim)
!  -K_ba^jr x L_rc^ik
       call DGEMM_('N','N',nadim,cdim,N,-one,ka(1,1,j),nadim,lc(1,ik),N,&
     &                     one,t3a,nadim)
!  -K_ca^ir x L_rb^jk   ! in the matrix as b,c,a
       call DGEMM_('T','T',bdim,ncdim,N,-one,lb(1,jk),N,kc(1,1,i),ncdim,&
     &                     zero,t3b,bdim)
!  -K_ca^kr x L_rb^ij
       call DGEMM_('T','T',bdim,ncdim,N,-one,lb(1,ij),N,kc(1,1,k),ncdim,&
     &                     one ,t3b,bdim)
!  K_ca^jr x L_rb^ik
       call DGEMM_('T','T',bdim,ncdim,N,one,lb(1,ik),N,kc(1,1,j),ncdim, &
     &                     one ,t3b,bdim)
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
!  K_cb^ir x L_ra^jk   ! in the matrix as b,c,a
       call DGEMM_('N','N',nbdim,adim,N,one,kb(1,1,i),nbdim,la(1,jk),N, &
     &                                  zero,t3a,nbdim)
!  K_cb^kr x L_ra^ij
       call DGEMM_('N','N',nbdim,adim,N,one,kb(1,1,k),nbdim,la(1,ij),N, &
     &                                  one ,t3a,nbdim)
!  -K_cb^jr x L_ra^ik
       call DGEMM_('N','N',nbdim,adim,N,-one,kb(1,1,j),nbdim,la(1,ik),N,&
     &                                  one ,t3a,nbdim)
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
! t3a  bac
! blok1   voa ka
         call DGEMM_('N','N',1,cdim,nadim, one,voa(1,ij),1,             &
     &  t3a,nadim,one,t1ac(k,1),noab)
         call DGEMM_('N','N',1,cdim,nadim, one,voa(1,jk),1,             &
     &  t3a,nadim,one,t1ac(i,1),noab)
         call DGEMM_('N','N',1,cdim,nadim,-one,voa(1,ik),1,             &
     &  t3a,nadim,one,t1ac(j,1),noab)
! blok 2 t3b bca
            call DGEMM_('N','T',1,bdim,ncdim,-one, voc(1,ij),1,         &
     &  t3b,bdim,one,t1ab(k,1),noab)
            call DGEMM_('N','T',1,bdim,ncdim,-one, voc(1,jk),1,         &
     &  t3b,bdim,one,t1ab(i,1),noab)
            call DGEMM_('N','T',1,bdim,ncdim,one,voc(1,ik),1,           &
     &  t3b,bdim,one,t1ab(j,1),noab)
            if(ifvo)then
             call DGEMM_('N','N',1,cdim,nadim, one,ka(1,i,j),1,         &
     &  t3a,nadim,one,t1bc(k,1),noab)
             call DGEMM_('N','N',1,cdim,nadim, one,ka(1,j,k),1,         &
     &  t3a,nadim,one,t1bc(i,1),noab)
             call DGEMM_('N','N',1,cdim,nadim,-one,ka(1,i,k),1,         &
     &  t3a,nadim,one,t1bc(j,1),noab)
            call DGEMM_('N','T',1,bdim,ncdim,-one, kc(1,i,j),1,         &
     &  t3b,bdim,one,t1bb(k,1),noab)
            call DGEMM_('N','T',1,bdim,ncdim,-one, kc(1,j,k),1,         &
     &  t3b,bdim,one,t1bb(i,1),noab)
            call DGEMM_('N','T',1,bdim,ncdim, one,kc(1,i,k),1,          &
     &  t3b,bdim,one,t1bb(j,1),noab)
             endif
! part 3 acb in t3b
            call transm(t3a,t3b,bdim,ncdim)
              call DGEMM_('N','T',1,adim,nbdim,one, vob(1,ij),1,        &
     &  t3b,adim,one,t1aa(k,1),noab)
               call DGEMM_('N','T',1,adim,nbdim,one, vob(1,jk),1,       &
     &  t3b,adim,one,t1aa(i,1),noab)
               call DGEMM_('N','T',1,adim,nbdim,-one,vob(1,ik),1,       &
     &  t3b,adim,one,t1aa(j,1),noab)
               if(ifvo)then
               call DGEMM_('N','T',1,adim,nbdim,one,  kb(1,i,j),1,      &
     &  t3b,adim,one,t1ba(k,1),noab)
               call DGEMM_('N','T',1,adim,nbdim,one,  kb(1,j,k),1,      &
     &  t3b,adim,one,t1ba(i,1),noab)
               call DGEMM_('N','T',1,adim,nbdim,-one, kb(1,i,k),1,      &
     &  t3b,adim,one,t1ba(j,1),noab)
               endif
         enddo !k
         enddo !j
         enddo !i
      return
! Avoid unused argument warnings
      if (.false.) then
        call Unused_real_array(mi)
        call Unused_real_array(mij)
        call Unused_integer(nuab)
      end if
      end
