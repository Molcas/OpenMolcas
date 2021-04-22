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
      subroutine t3_bt_aaa(nug,ka,kb,kc,la,mi,mij,adim,N,noab,
     $nuab,nnoab,lu,iasblock,nga,oeh,oep,enx,voa,t1a,t1b,t3a,t3b,ifvo)
      implicit none
      real*8 one,zero
      parameter (one=1.d0,zero=0.d0)
      integer a,b,c,aa,ab,bc,ac
      real*8 XX,YY,enx,den,dena,denb,denc
      integer ndim,adim,noab,nuab,i,j,k,nga,iasblock(3),lu(2),N
      integer ias,nga_offset,nug_offset,jk,ij,ik,nug,nnoab
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
         jk=0
         do j=2,i-1
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
