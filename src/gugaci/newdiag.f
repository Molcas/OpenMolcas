************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2008, Bingbing Suo                                     *
************************************************************************
C 04 Jul 2008 -BSUO- generalized davidson diagnolization routines
      subroutine gdavdiag(istate,mtidx,veistate)
#include "drt_h.fh"
      parameter (maxdimlu=1000,maxdimgit=10000)

c compute (h00-rouk)-1, and save the low triagular part of the
c inverse matrix at value_tmp(1:lent1-1)
      call cimtinverse(ndimt,veistate)

c h00 (h00-rouk)-1 * vector1(1:ndim)

      lent1=maxdimlu*maxdimlu
      lent2=lent1+ndimt
      value_lpext(lent1+1:lent2)=vector1(mtidx+1:mtidx+ndimt)
c subroutine dsymv in blas is used
c      call dsymv('l',ndimt,1.0,ndimt,value_lpext(1),1,
c     *           vector1(mtidx+1),1)
      call matmultv(value_lpext,ndimt,maxdimlu,
     *              value_lpext(lent1+1:lent1+ndimt),
     *              vector1(mtidx+1:ndimt))
c   from h00+1 to nci_dim
      do l=ndimt+1,nci_dim
         depff=vector2(l)-veistate
c               if(abs(depff).le.depc) depff=depc
         vector1(mtidx+l)=vector1(mtidx+l)/depff
      enddo

c Avoid unused argument warnings
      if (.false.) call Unused_integer(istate)
      end

      subroutine matmultv(a,n,np,x,y)
c matrix a(n,n), vector x(n),y(n)
c ax=y
      implicit real*8 (a-h,o-z)
      dimension a(np,np),x(n),y(n)

      y(1:n)=0.d0
      do i=1,n
        do j=1,n
          y(i)=y(i)+a(i,j)*x(j)
        enddo
      enddo

      return
c...end of matmulv
      end

      subroutine cimtinverse(ndimt,vrouk)
#include "drt_h.fh"
      common /scratch/ tmpdir,len_str
      character*256 tmpdir
!      character*256 filename

      parameter (maxdimlu=1000,maxdimgit=10000)

c      is=1
c      if(is.eq.1) then
c v'-v space is chosen as h0 part in gdvdd
c        ndimt=iw_sta(2,1)
c      else
c v'+d'+t'+s'-v space
c        ndimt=iw_sta(1,2)
c      endif

c      if(ndimt.le.maxdimlu) then
!        filename=tmpdir(1:len_str)//"/fort.23"
        lent=len_str+8
!        write(6,*) "new diag"
!        stop 1999
!        open(nf23,file=filename(1:lent),form="unformatted")
!        read(nf23) ndimt,lent,value_lpext(1:lent)
!        close(nf23)
        neh0=maxdimlu*maxdimlu
        num=2*neh0+lent
        if(num.gt.max_tmpvalue) then
          write(6,*) "no enough scracht vectors, error 1000"
     *               ,num,max_tmpvalue,neh0,lent,ndimt
#ifdef MOLPRO
#else
      call abend()
#endif
#ifdef _XIANEST_
#endif
!         call abend
!          stop 1000
        endif
        value_lpext(2*neh0+1:2*neh0+lent)=value_lpext(1:lent)
        call matinverse(value_lpext,
     *                  value_lpext(neh0+1:neh0+maxdimlu**2),ndimt,
     *                  maxdimlu,value_lpext(2*neh0+1:2*neh0+lent),lent,
     *                  vrouk)
c      endif

c      stop 1000
c the inverse matrix a is also a symmetric matrix,
c so we only need to keep the low triagular part of a.
c      call savelowtra(value_lpext(1),value_lpext(neh0+1),ndimt,
c     *                maxdimlu,neh0)
      value_lpext(1:neh0)=value_lpext(neh0+1:2*neh0)
      value_lpext(neh0+1:2*neh0)=0.d0
      return
      end

      subroutine matinverse(a,y,n,np,valvec,lent,vrouk)
      implicit real*8 (a-h,o-z)
      parameter (maxdimlu=1000,maxdimgit=10000)
      integer np, indx(n)
      dimension a(np,np),y(np,np),valvec(lent)
      parameter (zero=0.d0)
*i64  parameter (zero=0.e0)

      k=0
      do i=1,n
       do j=1,i-1
          k=k+1
          a(i,j)=valvec(k)
          a(j,i)=valvec(k)
        enddo
        k=k+1
        a(i,i)=valvec(k)-vrouk
      enddo
c      write(6,*) "ndim0 v-v,",np
c      do i=1,10
c        write(6,"(10(f12.6,1x))") a(i,1:10)
c      enddo
c      stop 888

!      temporay desiable next five line because these code is not used now
c      info = 0
c subroutine dgetrf in lapack is used
c      call dgetrf_lapack( n, n, a, n, indx, info )
c      info = 0
c      lwork = n*n
c subroutine dgetri in lapack is used
c      call dgetri_lapack( n, a, n, indx, y, lwork, info )

c      !y(1:np,1:np)=zero
c      !do i=1,n
c      !  y(i,i)=1.0
c      !enddo
c      !call ludcmp(a,n,np,indx,d)
c      !do i=1,n
c      !   call lubksb(a,n,np,indx,y(1,i))
c      !enddo

      return
      if (.false.) then
        call Unused_integer_array(indx)
        call Unused_real_array(y)
      endif
      end

      subroutine savelowtra(varry,a,ndimt,maxdimlu,neh0)
      implicit real*8 (a-h,o-z)
      dimension varry(neh0),a(maxdimlu,maxdimlu)
      parameter (dzero=0.d0)
*i64  parameter (dzero=0.e0)

      varry(1:neh0)=dzero
      k=0
      do i=1,ndimt
        do j=1,i-1
          k=k+1
          varry(k)=a(i,j)
        enddo
        k=k+1
        varry(k)=a(i,i)
      enddo

      return
c...end of savelowtra
      end
