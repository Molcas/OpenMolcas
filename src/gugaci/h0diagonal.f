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
      subroutine hymat_2(maxroot,minspace,ndim,kval,mroot,dcrita,eval,
     *                   vcm,indx,th,nxh,vb1,vb2,nxb,vad)
c     *****************************************************************
c     the sub. based on the algorithm of davidson (e.r. davison, j.
c     comp. phys. 17,87 (1975)) for searching the m-th eigenvalue and
c     eigenvector of the symmetric matrix . this sub. is an improved
c     version of the one given by j. weber , r. lacroix and g. wanner
c     in computers & chemistry , 4 , 55 , 1980 .
c     *****************************************************************
      implicit real*8 (a,c,d,e,f,p,h,w,v,t,r)
      !parameter (maxroot=10,minspace=40)
!      common /file_descript/nf1,nf2,nf3,nf4,nf7,nf8,nf9,nf10,
!     *                      nf11,nf13,nf15,nf20
      dimension indx(minspace),vd(minspace),ve(minspace),
     :   vu(minspace,minspace),vp(minspace*(minspace+1)/2),
     :   th(nxh),vb1(nxb),vb2(nxb),vad(ndim)
      dimension eval(maxroot),vcm(ndim*mroot),eeval(maxroot)
      dimension residvb(maxroot),valpha(maxroot),deff(maxroot)
      dimension ecrita(maxroot)
      data depc/1.0d-7/
c**************************************************************
c
c      write(6,*) 'generate vector vb2 from matrix a and vector vb1'
c
c**************************************************************

      ! debugging
      !vb1(1:ndim)=0.d0; vb2(1:ndim)=0.d0
      !vb1(5)=1.d0
      !call abprod2(ndim,1,kval,th,nxh,vb1,vb2,nxb,vad)
      !do i=1,ndim
      !  write(6,"(i8,f18.8)") i,vb2(i)
      !enddo
      !stop 888

      ecrita=1.e-8
      iterat=1
      mrsta=1
      call abprod2(ndim,1,kval,th,nxh,vb1,vb2,nxb,vad)
      call matrmk2(ndim,1,kval,indx,vp,vb1,vb2,nxb)

c==============================================================
      write(6,*)
      l=0
      do k=1,kval
        write(6,1112) (vp(i),i=l+1,l+k)
        l=l+k
      enddo
      write(6,*)
c==============================================================
c
c     write(6,*) 'diagonalization of matrix p(j,j)'
c
c==============================================================
200   iterat=iterat+1
      if(iterat.eq.200) then
        write(6,*) " h0 space fail to converged"
        write(6,*) " program stop"
        call abend()
      endif
      eeval(1:mroot)=eval(1:mroot)
      call hotred(minspace,kval,vp,vd,ve,vu)
      call qlcm(minspace,kval,vd,ve,vu)

      valpha(mrsta:mroot)=0.d0
      do m=mrsta,mroot
        eval(m)=vd(m)
        do k=1,kval
          tm=abs(vu(k,m))
         if(valpha(m).lt.tm) valpha(m)=tm   !max(vu(k,m),k=1,kval)
        enddo
        valpha(m)=1-valpha(m)*valpha(m)
        if(valpha(m).gt.depc) valpha(m)=sqrt(valpha(m))
        deff(m)  =abs(eval(m)-eeval(m))
      enddo
c****************************************************************
c
c     construction of the new approximate eigenvector vcm(n)'
c
c****************************************************************
      ijm=indx(mrsta)
      vcm(ijm+1:ndim*mroot)=0.0d0
      do k=1,kval
        ijb=indx(k)
        do m=mrsta,mroot
          ijm=indx(m)
          vukm=vu(k,m)
          do i=1,ndim
            vcm(ijm+i)=vcm(ijm+i)+vukm*vb1(ijb+i)
          enddo
        enddo
      enddo

      !do i=1,ndim
      !    write(6,"(i8,f18.8)") i, vcm(i)
      !enddo


      write(6,1113) iterat,kval,
     :        (m,eval(m),deff(m),m=mrsta,mroot)
1113  format(2i3,10(2x,i2,f14.8,f14.8))

      nroot=mroot-mrsta+1

      if(kval.eq.mroot*2) goto 10

      mrsta0=mrsta
      do m=mrsta0,mroot
c        if(deff(m).lt.ecrita(m).or.valpha(m).lt.dcrita) then     ! conv
        if(m.eq.mrsta.and.deff(m).lt.ecrita(m)) then              ! conv
c         if(valpha(m).lt.dcrita) then                           ! conve
         mrsta=mrsta+1
        endif
      enddo
c      mrsta=mrsta+mrooted
      nroot=mroot-mrsta+1
      if(mrsta.gt.mroot) then
        write(6,*)
        write(6,*)mroot,' roots are convegenced,after',iterat,' iterat'
        goto 300
      endif

10    mmspace=min(mroot*3+10,ndim)
!      nroot=1             ! bbs debug error?
      if(kval+nroot.gt.mmspace) then

c===== start  reset_basis ======================================

        do m=mrsta,kval
          indxm=indx(m)
          vb1(indxm+1:indxm+ndim)=0.0d0
        enddo
        do m=mrsta,mroot
          ijm=indx(m)
          do k=1,kval
            ijb=indx(k)
            vukm=vu(k,m)
            do l=1,ndim
              vb1(ijm+l)=vb1(ijm+l)+vukm*vb2(ijb+l)  ! h*cm-->vb1
            enddo
          enddo
        enddo
        ijm=indx(mrsta)
        vb2(ijm+1:ijm+ndim*nroot)=vb1(ijm+1:ijm+ndim*nroot)  ! h*cm-->vb
        vb1(ijm+1:ijm+ndim*nroot)=vcm(ijm+1:ijm+ndim*nroot)  !   cm-->vb
        residvb(mrsta:mroot)=0.d0

        mval=mroot
        do m=mrsta,mroot
          ijm=indx(m)
          mval=mval+1
          ijmb1=indx(mval)
          do l=1,ndim
            depcc= eval(m)-vad(l)
            if(abs(depcc).lt.depc) depcc=depc
            vb1(ijmb1+l)=(vb2(ijm+l)-vcm(ijm+l)*eval(m))/depcc
            residvb(m)=residvb(m)+vb1(ijmb1+l)*vb1(ijmb1+l)
          enddo
          call orthnor(ndim,mval,dcrita,vb1,nxb)
        enddo

        vb2(indx(mroot+1)+1:indx(kval+1)) =0.d0
        kval=mroot+nroot
        mval=mroot+1
        mn=mroot*(mroot+1)/2
        vp(1:mn)=0.0d0
        mn=0
        do m=1,mroot
          mn=mn+m
          vp(mn)=eval(m)
        enddo
        goto 100
c===== end  reset_basis ======================================
      endif

c
c     form the (j+1)-th approximate vector , vb1(n,j+1)
c
      jib1=indx(kval)
      jicm=indx(mrsta)

      residvb(mrsta:mroot)=0.d0
      do m=mrsta,mroot
        jib1=jib1+ndim
        do k=1,kval
          jib2=indx(k)
          do l=1,ndim
            vb1(jib1+l)=vb1(jib1+l)+vu(k,m)*vb2(jib2+l)
          enddo
        enddo

        do l=1,ndim
          depcc= eval(m)-vad(l)
          if(abs(depcc).lt.depc) depcc=depc
          vb1(jib1+l)=(vb1(jib1+l)-vcm(jicm+l)*eval(m))/depcc
          residvb(m)=residvb(m)+vb1(jib1+l)*vb1(jib1+l)
        enddo
        jicm=jicm+ndim
      enddo
      mval=kval+1
      do m=mrsta,mroot
        kval=kval+1
        call orthnor(ndim,kval,dcrita,vb1,nxb)
      enddo

100   call abprod2(ndim,mval,kval,th,nxh,vb1,vb2,nxb,vad)
      call matrmk2(ndim,mval,kval,indx,vp,vb1,vb2,nxb)

c=====  write out p_matrix =======================================
!      write(6,*)
!      write(nf2,*)
!      l=0
!      do k=1,kval
!        write(6 ,1112) (vp(i),i=l+1,l+k)
!        write(nf2,1112) (vp(i),i=l+1,l+k)
!        l=l+k
!     enddo
!      write(6,*)
c=====  write out p_matrix =======================================

      goto 200
c
300   continue

      ! copy ci vector to VB1
      do m=1,mroot
        ijm=indx(m)
        vb1(ijm+1:ijm+ndim)=vcm(ijm+1:ijm+ndim)
!       write(nf1) eval(m),(vcm(ijm+i),i=1,ndim)
!        write(6,"(5(1x,f18.9),1x,i2)") (vcm(ijm+i),i=1,4),vcm(35)
      enddo

1112  format(2x,20f14.8)

      return
      end

      subroutine matrmk2(n,k1,k2,indx,p,vb1,vb2,nxb)
      implicit real*8  (a-h,o-z)
      dimension p(465),indx(30),vb1(nxb),vb2(nxb)
      do 200 i=k1,k2
      iijj=i*(i-1)/2
      ij=indx(i)
      do 201 j=1,i
      ji=indx(j)
      p(iijj+j)=0.0d0
c--------------------------------------------------------------
      do 202 l=1,n
c-----------------------------------------------------------------
      p(iijj+j)=p(iijj+j)+vb1(ij+l)*vb2(ji+l)
202   continue
201   continue
200   continue
      return
      end

      subroutine abprod2(n,k1,k2,th,nxh,vb1,vb2,nxb,vad)
#include "drt_h.fh"
      dimension th(nxh),vb1(nxb),vb2(nxb),vad(n)
      !real*8,allocatable :: buff(:)
      !allocate(buff(n))
      ij=0
      do j=k1,k2
        ij=indx(j)
        do i=1,n
          vb2(ij+i)=vad(i)*vb1(ij+i)
        enddo
      enddo
c-------------------------------------------------------------------
      do 200 i=2,n
      mn=i*(i-1)/2
      do 201 j=k1,k2
      ij=indx(j)
      do 202 l=1,i-1
      vb2(ij+i)=vb2(ij+i)+th(mn+l)*vb1(ij+l)
      vb2(ij+l)=vb2(ij+l)+th(mn+l)*vb1(ij+i)
202   continue
201   continue
200   continue
c-------------------------------------------------------------------
      return
      end
c
      subroutine orthnor(n,j,dcrita,vb1,nxb)
#include "drt_h.fh"
      dimension vb1(nxb)
      ji=indx(j)
      if(j.eq.1) goto 150
      jm=j-1
      smax2=1.d10
120   smax1=0.0d0
      do 140 l=1,jm
      s=0.0d0
      ij=indx(l)
      do 130 i=1,n
      s=s+vb1(ij+i)*vb1(ji+i)
130   continue
      smax1=max(smax1,abs(s))
      do 141 i=1,n
        vb1(ji+i)=vb1(ji+i)-s*vb1(ij+i)
141   continue
140   continue

      if(smax1.lt.dcrita) goto 150
      if(smax1.gt.smax2) then
         write(6,*) 'dgnalization procedure is non-convergent.'
#ifdef MOLPRO
#else
      call abend()
#endif
#ifdef _XIANEST_
#endif
!        call abend
!         stop
      endif
      smax2=smax1
      goto 120
c     normalization of j-th eigenvector.
150   s=0.0d0
      do 160 i=1,n
      s=s+vb1(ji+i)*vb1(ji+i)
160   continue
      s=sqrt(s)
      do 170 i=1,n
      vb1(ji+i)=vb1(ji+i)/s
170   continue
      return
      end
c

      subroutine norm_a(n,av)  !bv:basis, av:vector for orth and norm
      real*8 av(n),s,ddot_,dcrita
      dcrita=1.0d-10
c     normalization of av_eigenvector.
      s=0.0d0
      s=ddot_(n,av,1,av,1)
      s=sqrt(s)
      s=max(s,dcrita)
      do i=1,n
        av(i)=av(i)/s
      enddo

      return
      end

      subroutine orth_ab(n,av,bv)  !bv:basis, av:vector for orth
      real*8 av(n),bv(n),s,ddot_
c     orthogonalization av,bv
      s=ddot_(n,av,1,bv,1)

      do i=1,n
        av(i)=av(i)-s*bv(i)
      enddo
      return
      end
