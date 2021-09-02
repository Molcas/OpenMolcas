!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2008, Bingbing Suo                                     *
!***********************************************************************
! 04 Jul 2008 -BSUO- generalized davidson diagnolization routines
      subroutine gdavdiag(istate,mtidx,veistate)
#include "drt_h.fh"
      parameter (maxdimlu=1000,maxdimgit=10000)

! compute (h00-rouk)-1, and save the low triagular part of the
! inverse matrix at value_tmp(1:lent1-1)
      call cimtinverse(ndimt,veistate)

! h00 (h00-rouk)-1 * vector1(1:ndim)

      lent1=maxdimlu*maxdimlu
      lent2=lent1+ndimt
      value_lpext(lent1+1:lent2)=vector1(mtidx+1:mtidx+ndimt)
! subroutine dsymv in blas is used
!      call dsymv('l',ndimt,1.0,ndimt,value_lpext(1),1,
!     *           vector1(mtidx+1),1)
      call matmultv(value_lpext,ndimt,maxdimlu,                         &
     &              value_lpext(lent1+1:lent1+ndimt),                   &
     &              vector1(mtidx+1:ndimt))
!   from h00+1 to nci_dim
      do l=ndimt+1,nci_dim
         depff=vector2(l)-veistate
!               if(abs(depff).le.depc) depff=depc
         vector1(mtidx+l)=vector1(mtidx+l)/depff
      enddo

! Avoid unused argument warnings
      if (.false.) call Unused_integer(istate)
      end

      subroutine matmultv(a,n,np,x,y)
! matrix a(n,n), vector x(n),y(n)
! ax=y
      implicit real*8 (a-h,o-z)
      dimension a(np,np),x(n),y(n)

      y(1:n)=0.d0
      do i=1,n
        do j=1,n
          y(i)=y(i)+a(i,j)*x(j)
        enddo
      enddo

      return
!...end of matmulv
      end

      subroutine cimtinverse(ndimt,vrouk)
#include "drt_h.fh"
#include "scratch.fh"
!      character*256 filename

      parameter (maxdimlu=1000,maxdimgit=10000)

!      is=1
!      if(is.eq.1) then
! v'-v space is chosen as h0 part in gdvdd
!        ndimt=iw_sta(2,1)
!      else
! v'+d'+t'+s'-v space
!        ndimt=iw_sta(1,2)
!      endif

!      if(ndimt.le.maxdimlu) then
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
          write(6,*) "no enough scracht vectors, error 1000"            &
     &               ,num,max_tmpvalue,neh0,lent,ndimt
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
        call matinverse(value_lpext,                                    &
     &                  value_lpext(neh0+1:neh0+maxdimlu**2),ndimt,     &
     &                  maxdimlu,value_lpext(2*neh0+1:2*neh0+lent),lent,&
     &                  vrouk)
!      endif

!      stop 1000
! the inverse matrix a is also a symmetric matrix,
! so we only need to keep the low triagular part of a.
!      call savelowtra(value_lpext(1),value_lpext(neh0+1),ndimt,
!     *                maxdimlu,neh0)
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
!i64  parameter (zero=0.e0)

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
!      write(6,*) "ndim0 v-v,",np
!      do i=1,10
!        write(6,"(10(f12.6,1x))") a(i,1:10)
!      enddo
!      stop 888

!      temporay desiable next five line because these code is not used now
!      info = 0
! subroutine dgetrf in lapack is used
!      call dgetrf_lapack( n, n, a, n, indx, info )
!      info = 0
!      lwork = n*n
! subroutine dgetri in lapack is used
!      call dgetri_lapack( n, a, n, indx, y, lwork, info )

!      !y(1:np,1:np)=zero
!      !do i=1,n
!      !  y(i,i)=1.0
!      !enddo
!      !call ludcmp(a,n,np,indx,d)
!      !do i=1,n
!      !   call lubksb(a,n,np,indx,y(1,i))
!      !enddo

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
!i64  parameter (dzero=0.e0)

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
!...end of savelowtra
      end
