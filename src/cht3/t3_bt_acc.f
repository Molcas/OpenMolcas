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
      subroutine t3_bt_acc(nug,ka,kb,kc,la,lc,mi,mij,adim,cdim,         &
     &N,noab,nuab,nnoab,lu,iasblock,nga,ngc,oeh,oepa,oepc,enx,voa,voc,  &
     &t1aa,t1ba,t1ac,t1bc,t3a,t3b,ifvo)
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
!
      sumt3=0.d0
      if(cdim.eq.1)return
      ncdim=cdim*(cdim-1)/2
      nadim=adim*cdim
      nug_offset=iasblock(1)*nug*(nug+1)/2
      ias=iasblock(2)*(nga-1)+1
      call multi_readir(la,nnoab*adim*N,lu(2),ias)
      ias=iasblock(2)*(ngc-1)+1
      call multi_readir(lc,nnoab*cdim*N,lu(2),ias)
! reads vvoo
      nga_offset=iasblock(3)*(nga*(nga-1)/2+ngc-1)+1
      ngc_offset=iasblock(3)*(ngc*(ngc-1)/2+ngc-1)+1
      ias=iasblock(2)*nug+nga_offset
      call multi_readir(voa,nnoab*nadim,lu(2),ias)
      ias=iasblock(2)*nug+ngc_offset
      call multi_readir(voc,nnoab*ncdim,lu(2),ias)
!
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
!  K_ab^ir x L_rc^jk   (stored as   c,b,a)
        call DGEMM_('T','T',cdim,nadim,N,one, lc(1,jk),N,               &
     &  ka(1,1,i),nadim,zero,t3a,cdim)
!  K_ab^kr x L_rc^ij
        call DGEMM_('T','T',cdim,nadim,N,one, lc(1,ij),N,ka(1,1,k),     &
     &              nadim,one,t3a,cdim)
!  -K_ab^jr x L_rc^ik
        call DGEMM_('T','T',cdim,nadim,N,-one,lc(1,ik),N,ka(1,1,j),     &
     &              nadim,one,t3a,cdim)
!
!  K_bc^ir x L_ra^jk
        call DGEMM_('N','N',ncdim,adim,N,one,kc(1,1,i),ncdim,la(1,jk),N,&
     &                     zero,t3b,ncdim)
!  K_bc^kr x L_ra^ij
        call DGEMM_('N','N',ncdim,adim,N,one,kc(1,1,k),ncdim,la(1,ij),N,&
     &                     one,t3b,ncdim)
!  -K_bc^jr x L_ra^ik
        call DGEMM_('N','N',ncdim,adim,N,-one,kc(1,1,j),ncdim,          &
     &        la(1,ik),N,                                               &
     &                     one,t3b,ncdim)
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
! t3a  bac
! blok1   voa ka
            call DGEMM_('N','T',1,cdim,nadim,-one,voa(1,ij),1,          &
     &  t3a,cdim,one,t1ac(k,1),noab)
            call DGEMM_('N','T',1,cdim,nadim,-one,voa(1,jk),1,          &
     &  t3a,cdim,one,t1ac(i,1),noab)
            call DGEMM_('N','T',1,cdim,nadim, one,voa(1,ik),1,          &
     &  t3a,cdim,one,t1ac(j,1),noab)
! blok 2 t3b bca
            call DGEMM_('N','N',1,adim,ncdim,one, voc(1,ij),1,          &
     &  t3b,ncdim,one,t1aa(k,1),noab)
            call DGEMM_('N','N',1,adim,ncdim,one, voc(1,jk),1,          &
     &  t3b,ncdim,one,t1aa(i,1),noab)
            call DGEMM_('N','N',1,adim,ncdim,-one,voc(1,ik),1,          &
     &  t3b,ncdim,one,t1aa(j,1),noab)
            if(ifvo)then
            call DGEMM_('N','T',1,cdim,nadim,-one,ka(1,i,j),1,          &
     &  t3a,cdim,one,t1bc(k,1),noab)
            call DGEMM_('N','T',1,cdim,nadim,-one,ka(1,j,k),1,          &
     &  t3a,cdim,one,t1bc(i,1),noab)
            call DGEMM_('N','T',1,cdim,nadim, one,ka(1,i,k),1,          &
     &  t3a,cdim,one,t1bc(j,1),noab)
            call DGEMM_('N','N',1,adim,ncdim,one, kc(1,i,j),1,          &
     &  t3b,ncdim,one,t1ba(k,1),noab)
            call DGEMM_('N','N',1,adim,ncdim,one, kc(1,j,k),1,          &
     &  t3b,ncdim,one,t1ba(i,1),noab)
            call DGEMM_('N','N',1,adim,ncdim,-one,kc(1,i,k),1,          &
     &  t3b,ncdim,one,t1ba(j,1),noab)
            endif
         enddo !k
         enddo !j
         enddo !i
      return
! Avoid unused argument warnings
      if (.false.) then
        call Unused_real_array(kb)
        call Unused_real_array(mi)
        call Unused_real_array(mij)
        call Unused_integer(nuab)
      end if
      end
