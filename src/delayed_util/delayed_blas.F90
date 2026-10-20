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
use link_blas, only: lb_caxpy
use Definitions, only: BLASInt, BLASR4
implicit none
complex(kind=BLASR4) :: ca
integer(kind=BLASInt) :: incx, incy, n
complex(kind=BLASR4) :: cx(*), cy(*)
call lb_caxpy(n,ca,cx,incx,cy,incy)
end subroutine caxpy

subroutine cgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
use link_blas, only: lb_cgemm
use Definitions, only: BLASInt, BLASR4
implicit none
complex(kind=BLASR4) :: alpha, beta
integer(kind=BLASInt) :: k, lda, ldb, ldc, m, n
character :: transa, transb
complex(kind=BLASR4) :: a(lda,*), b(ldb,*), c(ldc,*)
call lb_cgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
end subroutine cgemm

function dasum(n,dx,incx)
use link_blas, only: lb_dasum
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: incx, n
real(kind=BLASR8) :: dx(*)
real(kind=BLASR8) :: dasum
dasum = lb_dasum(n,dx,incx)
end function dasum

subroutine daxpy(n,da,dx,incx,dy,incy)
use link_blas, only: lb_daxpy
use Definitions, only: BLASInt, BLASR8
implicit none
real(kind=BLASR8) :: da
integer(kind=BLASInt) :: incx, incy, n
real(kind=BLASR8) :: dx(*), dy(*)
call lb_daxpy(n,da,dx,incx,dy,incy)
end subroutine daxpy

function dcabs1(z)
use link_blas, only: lb_dcabs1
use Definitions, only: BLASR8
implicit none
complex(kind=BLASR8) :: z
real(kind=BLASR8) :: dcabs1
dcabs1 = lb_dcabs1(z)
end function dcabs1

subroutine dcopy(n,dx,incx,dy,incy)
use link_blas, only: lb_dcopy
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: incx, incy, n
real(kind=BLASR8) :: dx(*), dy(*)
call lb_dcopy(n,dx,incx,dy,incy)
end subroutine dcopy

function ddot(n,x,dx,y,dy)
use link_blas, only: lb_ddot
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: n, dx, dy
real(kind=BLASR8) :: x(*), y(*)
real(kind=BLASR8) :: ddot
ddot = lb_ddot(n,x,dx,y,dy)
end function ddot

subroutine dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
use link_blas, only: lb_dgemm
use Definitions, only: BLASInt, BLASR8
implicit none
real(kind=BLASR8) :: alpha, beta
integer(kind=BLASInt) :: k, lda, ldb, ldc, m, n
character :: transa, transb
real(kind=BLASR8) :: a(lda,*), b(ldb,*), c(ldc,*)
call lb_dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
end subroutine dgemm

subroutine dgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
use link_blas, only: lb_dgemv
use Definitions, only: BLASInt, BLASR8
implicit none
real(kind=BLASR8) :: alpha, beta
integer(kind=BLASInt) :: incx, incy, lda, m, n
character :: trans
real(kind=BLASR8) :: a(lda,*), x(*), y(*)
call lb_dgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
end subroutine dgemv

subroutine dger(m,n,alpha,x,incx,y,incy,a,lda)
use link_blas, only: lb_dger
use Definitions, only: BLASInt, BLASR8
implicit none
real(kind=BLASR8) :: alpha
integer(kind=BLASInt) :: incx, incy, lda, m, n
real(kind=BLASR8) :: a(lda,*), x(*), y(*)
call lb_dger(m,n,alpha,x,incx,y,incy,a,lda)
end subroutine dger

function dnrm2(n,x,incx)
use link_blas, only: lb_dnrm2
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: incx, n
real(kind=BLASR8) :: x(*)
real(kind=BLASR8) :: dnrm2
dnrm2 = lb_dnrm2(n,x,incx)
end function dnrm2

subroutine drot(n,dx,incx,dy,incy,c,s)
use link_blas, only: lb_drot
use Definitions, only: BLASInt, BLASR8
implicit none
real(kind=BLASR8) :: c, s
integer(kind=BLASInt) :: incx, incy, n
real(kind=BLASR8) :: dx(*), dy(*)
call lb_drot(n,dx,incx,dy,incy,c,s)
end subroutine drot

subroutine dscal(n,da,dx,incx)
use link_blas, only: lb_dscal
use Definitions, only: BLASInt, BLASR8
implicit none
real(kind=BLASR8) :: da
integer(kind=BLASInt) :: incx, n
real(kind=BLASR8) :: dx(*)
call lb_dscal(n,da,dx,incx)
end subroutine dscal

subroutine dspmv(uplo,n,alpha,ap,x,incx,beta,y,incy)
use link_blas, only: lb_dspmv
use Definitions, only: BLASInt, BLASR8
implicit none
real(kind=BLASR8) :: alpha, beta
integer(kind=BLASInt) :: incx, incy, n
character :: uplo
real(kind=BLASR8) :: ap(*), x(*), y(*)
call lb_dspmv(uplo,n,alpha,ap,x,incx,beta,y,incy)
end subroutine dspmv

subroutine dspr(uplo,n,alpha,x,incx,ap)
use link_blas, only: lb_dspr
use Definitions, only: BLASInt, BLASR8
implicit none
real(kind=BLASR8) :: alpha
integer(kind=BLASInt) :: incx, n
character :: uplo
real(kind=BLASR8) :: ap(*), x(*)
call lb_dspr(uplo,n,alpha,x,incx,ap)
end subroutine dspr

subroutine dspr2(uplo,n,alpha,x,incx,y,incy,ap)
use link_blas, only: lb_dspr2
use Definitions, only: BLASInt, BLASR8
implicit none
real(kind=BLASR8) :: alpha
integer(kind=BLASInt) :: incx, incy, n
character :: uplo
real(kind=BLASR8) :: ap(*), x(*), y(*)
call lb_dspr2(uplo,n,alpha,x,incx,y,incy,ap)
end subroutine dspr2

subroutine dswap(n,dx,incx,dy,incy)
use link_blas, only: lb_dswap
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: incx, incy, n
real(kind=BLASR8) :: dx(*), dy(*)
call lb_dswap(n,dx,incx,dy,incy)
end subroutine dswap

subroutine dsymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
use link_blas, only: lb_dsymm
use Definitions, only: BLASInt, BLASR8
implicit none
real(kind=BLASR8) :: alpha, beta
integer(kind=BLASInt) :: lda, ldb, ldc, m, n
character :: side, uplo
real(kind=BLASR8) :: a(lda,*), b(ldb,*), c(ldc,*)
call lb_dsymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
end subroutine dsymm

subroutine dsymv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
use link_blas, only: lb_dsymv
use Definitions, only: BLASInt, BLASR8
implicit none
real(kind=BLASR8) :: alpha, beta
integer(kind=BLASInt) :: incx, incy, lda, n
character :: uplo
real(kind=BLASR8) :: a(lda,*), x(*), y(*)
call lb_dsymv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
end subroutine dsymv

subroutine dsyr2(uplo,n,alpha,x,incx,y,incy,a,lda)
use link_blas, only: lb_dsyr2
use Definitions, only: BLASInt, BLASR8
implicit none
real(kind=BLASR8) :: alpha
integer(kind=BLASInt) :: incx, incy, lda, n
character :: uplo
real(kind=BLASR8) :: a(lda,*), x(*), y(*)
call lb_dsyr2(uplo,n,alpha,x,incx,y,incy,a,lda)
end subroutine dsyr2

subroutine dsyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
use link_blas, only: lb_dsyr2k
use Definitions, only: BLASInt, BLASR8
implicit none
real(kind=BLASR8) :: alpha, beta
integer(kind=BLASInt) :: k, lda, ldb, ldc, n
character :: trans, uplo
real(kind=BLASR8) :: a(lda,*), b(ldb,*), c(ldc,*)
call lb_dsyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
end subroutine dsyr2k

subroutine dsyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
use link_blas, only: lb_dsyrk
use Definitions, only: BLASInt, BLASR8
implicit none
real(kind=BLASR8) :: alpha, beta
integer(kind=BLASInt) :: k, lda, ldc, n
character :: trans, uplo
real(kind=BLASR8) :: a(lda,*), c(ldc,*)
call lb_dsyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
end subroutine dsyrk

subroutine dtpmv(uplo,trans,diag,n,ap,x,incx)
use link_blas, only: lb_dtpmv
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: incx, n
character :: diag, trans, uplo
real(kind=BLASR8) :: ap(*), x(*)
call lb_dtpmv(uplo,trans,diag,n,ap,x,incx)
end subroutine dtpmv

subroutine dtpsv(uplo,trans,diag,n,ap,x,incx)
use link_blas, only: lb_dtpsv
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: incx, n
character :: diag, trans, uplo
real(kind=BLASR8) :: ap(*), x(*)
call lb_dtpsv(uplo,trans,diag,n,ap,x,incx)
end subroutine dtpsv

subroutine dtrmm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
use link_blas, only: lb_dtrmm
use Definitions, only: BLASInt, BLASR8
implicit none
real(kind=BLASR8) :: alpha
integer(kind=BLASInt) :: lda, ldb, m, n
character :: diag, side, transa, uplo
real(kind=BLASR8) :: a(lda,*), b(ldb,*)
call lb_dtrmm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
end subroutine dtrmm

subroutine dtrmv(uplo,trans,diag,n,a,lda,x,incx)
use link_blas, only: lb_dtrmv
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: incx, lda, n
character :: diag, trans, uplo
real(kind=BLASR8) :: a(lda,*), x(*)
call lb_dtrmv(uplo,trans,diag,n,a,lda,x,incx)
end subroutine dtrmv

subroutine dtrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
use link_blas, only: lb_dtrsm
use Definitions, only: BLASInt, BLASR8
implicit none
real(kind=BLASR8) :: alpha
integer(kind=BLASInt) :: lda, ldb, m, n
character :: diag, side, transa, uplo
real(kind=BLASR8) :: a(lda,*), b(ldb,*)
call lb_dtrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
end subroutine dtrsm

subroutine dtrsv(uplo,trans,diag,n,a,lda,x,incx)
use link_blas, only: lb_dtrsv
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: incx, lda, n
character :: diag, trans, uplo
real(kind=BLASR8) :: a(lda,*), x(*)
call lb_dtrsv(uplo,trans,diag,n,a,lda,x,incx)
end subroutine dtrsv

function dznrm2(n,x,incx)
use link_blas, only: lb_dznrm2
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: incx, n
complex(kind=BLASR8) :: x(*)
real(kind=BLASR8) :: dznrm2
dznrm2 = lb_dznrm2(n,x,incx)
end function dznrm2

function idamax(n,dx,incx)
use link_blas, only: lb_idamax
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: incx, n
real(kind=BLASR8) :: dx(*)
integer(kind=BLASInt) :: idamax
idamax = lb_idamax(n,dx,incx)
end function idamax

function lsame(ca,cb)
use link_blas, only: lb_lsame
implicit none
character :: ca, cb
logical :: lsame
lsame = lb_lsame(ca,cb)
end function lsame

subroutine saxpy(n,sa,sx,incx,sy,incy)
use link_blas, only: lb_saxpy
use Definitions, only: BLASInt, BLASR4
implicit none
real(kind=BLASR4) :: sa
integer(kind=BLASInt) :: incx, incy, n
real(kind=BLASR4) :: sx(*), sy(*)
call lb_saxpy(n,sa,sx,incx,sy,incy)
end subroutine saxpy

function scabs1(z)
use link_blas, only: lb_scabs1
use Definitions, only: BLASR4
implicit none
complex(kind=BLASR4) :: z
real(kind=BLASR4) :: scabs1
scabs1 = lb_scabs1(z)
end function scabs1

subroutine scopy(n,sx,incx,sy,incy)
use link_blas, only: lb_scopy
use Definitions, only: BLASInt, BLASR4
implicit none
integer(kind=BLASInt) :: incx, incy, n
real(kind=BLASR4) :: sx(*), sy(*)
call lb_scopy(n,sx,incx,sy,incy)
end subroutine scopy

subroutine sgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
use link_blas, only: lb_sgemm
use Definitions, only: BLASInt, BLASR4
implicit none
real(kind=BLASR4) :: alpha, beta
integer(kind=BLASInt) :: k, lda, ldb, ldc, m, n
character :: transa, transb
real(kind=BLASR4) :: a(lda,*), b(ldb,*), c(ldc,*)
call lb_sgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
end subroutine sgemm

subroutine xerbla(srname,info)
use link_blas, only: lb_xerbla
use Definitions, only: BLASInt
implicit none
character(len=*) :: srname
integer(kind=BLASInt) :: info
call lb_xerbla(srname,info)
end subroutine xerbla

subroutine zaxpy(n,za,zx,incx,zy,incy)
use link_blas, only: lb_zaxpy
use Definitions, only: BLASInt, BLASR8
implicit none
complex(kind=BLASR8) :: za
integer(kind=BLASInt) :: incx, incy, n
complex(kind=BLASR8) :: zx(*), zy(*)
call lb_zaxpy(n,za,zx,incx,zy,incy)
end subroutine zaxpy

subroutine zcopy(n,zx,incx,zy,incy)
use link_blas, only: lb_zcopy
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: incx, incy, n
complex(kind=BLASR8) :: zx(*), zy(*)
call lb_zcopy(n,zx,incx,zy,incy)
end subroutine zcopy

function zdotc(n,zx,incx,zy,incy)
use link_blas, only: lb_zdotc
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: incx, incy, n
complex(kind=BLASR8) :: zx(*), zy(*)
complex(kind=BLASR8) :: zdotc
zdotc = lb_zdotc(n,zx,incx,zy,incy)
end function zdotc

subroutine zdscal(n,da,zx,incx)
use link_blas, only: lb_zdscal
use Definitions, only: BLASInt, BLASR8
implicit none
real(kind=BLASR8) :: da
integer(kind=BLASInt) :: incx, n
complex(kind=BLASR8) :: zx(*)
call lb_zdscal(n,da,zx,incx)
end subroutine zdscal

subroutine zdrot(n,cx,incx,cy,incy,c,s)
use link_blas, only: lb_zdrot
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: incx, incy, n
real(kind=BLASR8) :: c, s
complex(kind=BLASR8) :: cx(*), cy(*)
call lb_zdrot(n,cx,incx,cy,incy,c,s)
end subroutine zdrot

subroutine zgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
use link_blas, only: lb_zgemm
use Definitions, only: BLASInt, BLASR8
implicit none
complex(kind=BLASR8) :: alpha, beta
integer(kind=BLASInt) :: k, lda, ldb, ldc, m, n
character :: transa, transb
complex(kind=BLASR8) :: a(lda,*), b(ldb,*), c(ldc,*)
call lb_zgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
end subroutine zgemm

subroutine zgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
use link_blas, only: lb_zgemv
use Definitions, only: BLASInt, BLASR8
implicit none
complex(kind=BLASR8) :: alpha, beta
integer(kind=BLASInt) :: incx, incy, lda, m, n
character :: trans
complex(kind=BLASR8) :: a(lda,*), x(*), y(*)
call lb_zgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
end subroutine zgemv

subroutine zgerc(m,n,alpha,x,incx,y,incy,a,lda)
use link_blas, only: lb_zgerc
use Definitions, only: BLASInt, BLASR8
implicit none
complex(kind=BLASR8) :: alpha
integer(kind=BLASInt) :: incx, incy, lda, m, n
complex(kind=BLASR8) :: a(lda,*), x(*), y(*)
call lb_zgerc(m,n,alpha,x,incx,y,incy,a,lda)
end subroutine zgerc

subroutine zhemv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
use link_blas, only: lb_zhemv
use Definitions, only: BLASInt, BLASR8
implicit none
complex(kind=BLASR8) :: alpha, beta
integer(kind=BLASInt) :: incx, incy, lda, n
character :: uplo
complex(kind=BLASR8) :: a(lda,*), x(*), y(*)
call lb_zhemv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
end subroutine zhemv

subroutine zher2(uplo,n,alpha,x,incx,y,incy,a,lda)
use link_blas, only: lb_zher2
use Definitions, only: BLASInt, BLASR8
implicit none
complex(kind=BLASR8) :: alpha
integer(kind=BLASInt) :: incx, incy, lda, n
character :: uplo
complex(kind=BLASR8) :: a(lda,*), x(*), y(*)
call lb_zher2(uplo,n,alpha,x,incx,y,incy,a,lda)
end subroutine zher2

subroutine zher2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
use link_blas, only: lb_zher2k
use Definitions, only: BLASInt, BLASR8
implicit none
complex(kind=BLASR8) :: alpha
real(kind=BLASR8) :: beta
integer(kind=BLASInt) :: k, lda, ldb, ldc, n
character :: trans, uplo
complex(kind=BLASR8) :: a(lda,*), b(ldb,*), c(ldc,*)
call lb_zher2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
end subroutine zher2k

subroutine zhpmv(uplo,n,alpha,ap,x,incx,beta,y,incy)
use link_blas, only: lb_zhpmv
use Definitions, only: BLASInt, BLASR8
implicit none
complex(kind=BLASR8) :: alpha, beta
integer(kind=BLASInt) :: incx, incy, n
character :: uplo
complex(kind=BLASR8) :: ap(*), x(*), y(*)
call lb_zhpmv(uplo,n,alpha,ap,x,incx,beta,y,incy)
end subroutine zhpmv

subroutine zhpr2(uplo,n,alpha,x,incx,y,incy,ap)
use link_blas, only: lb_zhpr2
use Definitions, only: BLASInt, BLASR8
implicit none
complex(kind=BLASR8) :: alpha
integer(kind=BLASInt) :: incx, incy, n
character :: uplo
complex(kind=BLASR8) :: ap(*), x(*), y(*)
call lb_zhpr2(uplo,n,alpha,x,incx,y,incy,ap)
end subroutine zhpr2

subroutine zscal(n,za,zx,incx)
use link_blas, only: lb_zscal
use Definitions, only: BLASInt, BLASR8
implicit none
complex(kind=BLASR8) :: za
integer(kind=BLASInt) :: incx, n
complex(kind=BLASR8) :: zx(*)
call lb_zscal(n,za,zx,incx)
end subroutine zscal

subroutine zswap(n,zx,incx,zy,incy)
use link_blas, only: lb_zswap
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: incx, incy, n
complex(kind=BLASR8) :: zx(*), zy(*)
call lb_zswap(n,zx,incx,zy,incy)
end subroutine zswap

subroutine ztrmm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
use link_blas, only: lb_ztrmm
use Definitions, only: BLASInt, BLASR8
implicit none
complex(kind=BLASR8) :: alpha
integer(kind=BLASInt) :: lda, ldb, m, n
character :: diag, side, transa, uplo
complex(kind=BLASR8) :: a(lda,*), b(ldb,*)
call lb_ztrmm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
end subroutine ztrmm

subroutine ztrmv(uplo,trans,diag,n,a,lda,x,incx)
use link_blas, only: lb_ztrmv
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: incx, lda, n
character :: diag, trans, uplo
complex(kind=BLASR8) :: a(lda,*), x(*)
call lb_ztrmv(uplo,trans,diag,n,a,lda,x,incx)
end subroutine ztrmv

subroutine ztrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
use link_blas, only: lb_ztrsm
use Definitions, only: BLASInt, BLASR8
implicit none
character :: side, uplo, transa, diag
integer(kind=BLASInt) :: m, n, lda, ldb
complex(kind=BLASR8) :: alpha, a(lda,*), b(ldb,*)
call lb_ztrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
end subroutine ztrsm
