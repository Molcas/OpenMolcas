!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!**********************************************************************

subroutine caxpy(n,ca,cx,incx,cy,incy)
  use link_blas
  implicit none
  complex :: ca
  integer :: incx,incy,n
  complex :: cx(*),cy(*)
  call lb_caxpy(n,ca,cx,incx,cy,incy)
end subroutine caxpy

subroutine cgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
  use link_blas
  implicit none
  complex :: alpha,beta
  integer :: k,lda,ldb,ldc,m,n
  character :: transa,transb
  complex :: a(lda,*),b(ldb,*),c(ldc,*)
  call lb_cgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
end subroutine cgemm

function dasum(n,dx,incx)
  use link_blas
  implicit none
  integer :: incx,n
  real*8 :: dx(*)
  real*8 :: dasum
  dasum=lb_dasum(n,dx,incx)
end function dasum

subroutine daxpy(n,da,dx,incx,dy,incy)
  use link_blas
  implicit none
  real*8 :: da
  integer :: incx,incy,n
  real*8 :: dx(*),dy(*)
  call lb_daxpy(n,da,dx,incx,dy,incy)
end subroutine daxpy

function dcabs1(z)
  use link_blas
  implicit none
  complex*16 :: z
  real*8 :: dcabs1
  dcabs1=lb_dcabs1(z)
end function dcabs1

subroutine dcopy(n,dx,incx,dy,incy)
  use link_blas
  implicit none
  integer :: incx,incy,n
  real*8 :: dx(*),dy(*)
  call lb_dcopy(n,dx,incx,dy,incy)
end subroutine dcopy

function ddot(n,x,dx,y,dy)
  use link_blas
  implicit none
  integer :: n,dx,dy
  real*8 :: x(*),y(*)
  real*8 :: ddot
  ddot=lb_ddot(n,x,dx,y,dy)
end function ddot

subroutine dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
  use link_blas
  implicit none
  real*8 :: alpha,beta
  integer :: k,lda,ldb,ldc,m,n
  character :: transa,transb
  real*8 :: a(lda,*),b(ldb,*),c(ldc,*)
  call lb_dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
end subroutine dgemm

subroutine dgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
  use link_blas
  implicit none
  real*8 :: alpha,beta
  integer :: incx,incy,lda,m,n
  character :: trans
  real*8 :: a(lda,*),x(*),y(*)
  call lb_dgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
end subroutine dgemv

subroutine dger(m,n,alpha,x,incx,y,incy,a,lda)
  use link_blas
  implicit none
  real*8 :: alpha
  integer :: incx,incy,lda,m,n
  real*8 :: a(lda,*),x(*),y(*)
  call lb_dger(m,n,alpha,x,incx,y,incy,a,lda)
end subroutine dger

function dnrm2(n,x,incx)
  use link_blas
  implicit none
  integer :: incx,n
  real*8 :: x(*)
  real*8 :: dnrm2
  dnrm2=lb_dnrm2(n,x,incx)
end function dnrm2

subroutine drot(n,dx,incx,dy,incy,c,s)
  use link_blas
  implicit none
  real*8 :: c,s
  integer :: incx,incy,n
  real*8 :: dx(*),dy(*)
  call lb_drot(n,dx,incx,dy,incy,c,s)
end subroutine drot

subroutine dscal(n,da,dx,incx)
  use link_blas
  implicit none
  real*8 :: da
  integer :: incx,n
  real*8 :: dx(*)
  call lb_dscal(n,da,dx,incx)
end subroutine dscal

subroutine dspmv(uplo,n,alpha,ap,x,incx,beta,y,incy)
  use link_blas
  implicit none
  real*8 :: alpha,beta
  integer :: incx,incy,n
  character :: uplo
  real*8 :: ap(*),x(*),y(*)
  call lb_dspmv(uplo,n,alpha,ap,x,incx,beta,y,incy)
end subroutine dspmv

subroutine dspr(uplo,n,alpha,x,incx,ap)
  use link_blas
  implicit none
  real*8 :: alpha
  integer :: incx,n
  character :: uplo
  real*8 :: ap(*),x(*)
  call lb_dspr(uplo,n,alpha,x,incx,ap)
end subroutine dspr

subroutine dspr2(uplo,n,alpha,x,incx,y,incy,ap)
  use link_blas
  implicit none
  real*8 :: alpha
  integer :: incx,incy,n
  character :: uplo
  real*8 :: ap(*),x(*),y(*)
  call lb_dspr2(uplo,n,alpha,x,incx,y,incy,ap)
end subroutine dspr2

subroutine dswap(n,dx,incx,dy,incy)
  use link_blas
  implicit none
  integer :: incx,incy,n
  real*8 :: dx(*),dy(*)
  call lb_dswap(n,dx,incx,dy,incy)
end subroutine dswap

subroutine dsymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
  use link_blas
  implicit none
  real*8 :: alpha,beta
  integer :: lda,ldb,ldc,m,n
  character :: side,uplo
  real*8 :: a(lda,*),b(ldb,*),c(ldc,*)
  call lb_dsymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
end subroutine dsymm

subroutine dsymv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
  use link_blas
  implicit none
  real*8 :: alpha,beta
  integer :: incx,incy,lda,n
  character :: uplo
  real*8 :: a(lda,*),x(*),y(*)
  call lb_dsymv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
end subroutine dsymv

subroutine dsyr2(uplo,n,alpha,x,incx,y,incy,a,lda)
  use link_blas
  implicit none
  real*8 :: alpha
  integer :: incx,incy,lda,n
  character :: uplo
  real*8 :: a(lda,*),x(*),y(*)
  call lb_dsyr2(uplo,n,alpha,x,incx,y,incy,a,lda)
end subroutine dsyr2

subroutine dsyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
  use link_blas
  implicit none
  real*8 :: alpha,beta
  integer :: k,lda,ldb,ldc,n
  character :: trans,uplo
  real*8 :: a(lda,*),b(ldb,*),c(ldc,*)
  call lb_dsyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
end subroutine dsyr2k

subroutine dsyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
  use link_blas
  implicit none
  real*8 :: alpha,beta
  integer :: k,lda,ldc,n
  character :: trans,uplo
  real*8 :: a(lda,*),c(ldc,*)
  call lb_dsyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
end subroutine dsyrk

subroutine dtpmv(uplo,trans,diag,n,ap,x,incx)
  use link_blas
  implicit none
  integer :: incx,n
  character :: diag,trans,uplo
  real*8 :: ap(*),x(*)
  call lb_dtpmv(uplo,trans,diag,n,ap,x,incx)
end subroutine dtpmv

subroutine dtpsv(uplo,trans,diag,n,ap,x,incx)
  use link_blas
  implicit none
  integer :: incx,n
  character :: diag,trans,uplo
  real*8 :: ap(*),x(*)
  call lb_dtpsv(uplo,trans,diag,n,ap,x,incx)
end subroutine dtpsv

subroutine dtrmm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
  use link_blas
  implicit none
  real*8 :: alpha
  integer :: lda,ldb,m,n
  character :: diag,side,transa,uplo
  real*8 :: a(lda,*),b(ldb,*)
  call lb_dtrmm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
end subroutine dtrmm

subroutine dtrmv(uplo,trans,diag,n,a,lda,x,incx)
  use link_blas
  implicit none
  integer :: incx,lda,n
  character :: diag,trans,uplo
  real*8 :: a(lda,*),x(*)
  call lb_dtrmv(uplo,trans,diag,n,a,lda,x,incx)
end subroutine dtrmv

subroutine dtrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
  use link_blas
  implicit none
  real*8 :: alpha
  integer :: lda,ldb,m,n
  character :: diag,side,transa,uplo
  real*8 :: a(lda,*),b(ldb,*)
  call lb_dtrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
end subroutine dtrsm

subroutine dtrsv(uplo,trans,diag,n,a,lda,x,incx)
  use link_blas
  implicit none
  integer :: incx,lda,n
  character :: diag,trans,uplo
  real*8 :: a(lda,*),x(*)
  call lb_dtrsv(uplo,trans,diag,n,a,lda,x,incx)
end subroutine dtrsv

function dznrm2(n,x,incx)
  use link_blas
  implicit none
  integer :: incx,n
  complex*16 :: x(*)
  real*8 :: dznrm2
  dznrm2=lb_dznrm2(n,x,incx)
end function dznrm2

function idamax(n,dx,incx)
  use link_blas
  implicit none
  integer :: incx,n
  real*8 :: dx(*)
  integer :: idamax
  idamax=lb_idamax(n,dx,incx)
end function idamax

function lsame(ca,cb)
  use link_blas
  implicit none
  character :: ca,cb
  logical :: lsame
  lsame=lb_lsame(ca,cb)
end function lsame

subroutine saxpy(n,sa,sx,incx,sy,incy)
  use link_blas
  implicit none
  real :: sa
  integer :: incx,incy,n
  real :: sx(*),sy(*)
  call lb_saxpy(n,sa,sx,incx,sy,incy)
end subroutine saxpy

function scabs1(z)
  use link_blas
  implicit none
  complex :: z
  real :: scabs1
  scabs1=lb_scabs1(z)
end function scabs1

subroutine scopy(n,sx,incx,sy,incy)
  use link_blas
  implicit none
  integer :: incx,incy,n
  real :: sx(*),sy(*)
  call lb_scopy(n,sx,incx,sy,incy)
end subroutine scopy

subroutine sgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
  use link_blas
  implicit none
  real :: alpha,beta
  integer :: k,lda,ldb,ldc,m,n
  character :: transa,transb
  real :: a(lda,*),b(ldb,*),c(ldc,*)
  call lb_sgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
end subroutine sgemm

subroutine xerbla(srname,info)
  use link_blas
  implicit none
  character(len=*) :: srname
  integer :: info
  call lb_xerbla(srname,info)
end subroutine xerbla

subroutine zaxpy(n,za,zx,incx,zy,incy)
  use link_blas
  implicit none
  complex*16 :: za
  integer :: incx,incy,n
  complex*16 :: zx(*),zy(*)
  call lb_zaxpy(n,za,zx,incx,zy,incy)
end subroutine zaxpy

subroutine zcopy(n,zx,incx,zy,incy)
  use link_blas
  implicit none
  integer :: incx,incy,n
  complex*16 :: zx(*),zy(*)
  call lb_zcopy(n,zx,incx,zy,incy)
end subroutine zcopy

function zdotc(n,zx,incx,zy,incy)
  use link_blas
  implicit none
  integer :: incx,incy,n
  complex*16 :: zx(*),zy(*)
  complex*16 :: zdotc
  zdotc=lb_zdotc(n,zx,incx,zy,incy)
end function zdotc

subroutine zdscal(n,da,zx,incx)
  use link_blas
  implicit none
  real*8 :: da
  integer :: incx,n
  complex*16 :: zx(*)
  call lb_zdscal(n,da,zx,incx)
end subroutine zdscal

subroutine zgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
  use link_blas
  implicit none
  complex*16 :: alpha,beta
  integer :: k,lda,ldb,ldc,m,n
  character :: transa,transb
  complex*16 :: a(lda,*),b(ldb,*),c(ldc,*)
  call lb_zgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
end subroutine zgemm

subroutine zgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
  use link_blas
  implicit none
  complex*16 :: alpha,beta
  integer :: incx,incy,lda,m,n
  character :: trans
  complex*16 :: a(lda,*),x(*),y(*)
  call lb_zgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
end subroutine zgemv

subroutine zgerc(m,n,alpha,x,incx,y,incy,a,lda)
  use link_blas
  implicit none
  complex*16 :: alpha
  integer :: incx,incy,lda,m,n
  complex*16 :: a(lda,*),x(*),y(*)
  call lb_zgerc(m,n,alpha,x,incx,y,incy,a,lda)
end subroutine zgerc

subroutine zhemv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
  use link_blas
  implicit none
  complex*16 :: alpha,beta
  integer :: incx,incy,lda,n
  character :: uplo
  complex*16 :: a(lda,*),x(*),y(*)
  call lb_zhemv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
end subroutine zhemv

subroutine zher2(uplo,n,alpha,x,incx,y,incy,a,lda)
  use link_blas
  implicit none
  complex*16 :: alpha
  integer :: incx,incy,lda,n
  character :: uplo
  complex*16 :: a(lda,*),x(*),y(*)
  call lb_zher2(uplo,n,alpha,x,incx,y,incy,a,lda)
end subroutine zher2

subroutine zher2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
  use link_blas
  implicit none
  complex*16 :: alpha
  real*8 :: beta
  integer :: k,lda,ldb,ldc,n
  character :: trans,uplo
  complex*16 :: a(lda,*),b(ldb,*),c(ldc,*)
  call lb_zher2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
end subroutine zher2k

subroutine zhpmv(uplo,n,alpha,ap,x,incx,beta,y,incy)
  use link_blas
  implicit none
  complex*16 :: alpha,beta
  integer :: incx,incy,n
  character :: uplo
  complex*16 :: ap(*),x(*),y(*)
  call lb_zhpmv(uplo,n,alpha,ap,x,incx,beta,y,incy)
end subroutine zhpmv

subroutine zhpr2(uplo,n,alpha,x,incx,y,incy,ap)
  use link_blas
  implicit none
  complex*16 :: alpha
  integer :: incx,incy,n
  character :: uplo
  complex*16 :: ap(*),x(*),y(*)
  call lb_zhpr2(uplo,n,alpha,x,incx,y,incy,ap)
end subroutine zhpr2

subroutine zscal(n,za,zx,incx)
  use link_blas
  implicit none
  complex*16 :: za
  integer :: incx,n
  complex*16 :: zx(*)
  call lb_zscal(n,za,zx,incx)
end subroutine zscal

subroutine zswap(n,zx,incx,zy,incy)
  use link_blas
  implicit none
  integer :: incx,incy,n
  complex*16 :: zx(*),zy(*)
  call lb_zswap(n,zx,incx,zy,incy)
end subroutine zswap

subroutine ztrmm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
  use link_blas
  implicit none
  complex*16 :: alpha
  integer :: lda,ldb,m,n
  character :: diag,side,transa,uplo
  complex*16 :: a(lda,*),b(ldb,*)
  call lb_ztrmm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
end subroutine ztrmm

subroutine ztrmv(uplo,trans,diag,n,a,lda,x,incx)
  use link_blas
  implicit none
  integer :: incx,lda,n
  character :: diag,trans,uplo
  complex*16 :: a(lda,*),x(*)
  call lb_ztrmv(uplo,trans,diag,n,a,lda,x,incx)
end subroutine ztrmv

