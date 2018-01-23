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
! Molcas BLAS wrappers
! These wrappers handle lp64/ilp64 interfaces and various other
! requirements related to e.g. matrix size on GPUs. Any BLAS call in
! Molcas should use these wrappers.

! List of BLAS wrappers currently implemented:
! drot_
! dswap_
! dscal_
! dcopy_
! zcopy_
! scopy_
! daxpy_
! ddot_
! dnrm2_
! dasum_
! idamax_
! dgemm_
! dspmv_
! dgemv_
! dznrm2_
! zgemm_

! Set the appropriate integer size of the library interface:
#ifdef LINALG_I4
#  define BLASINT INTEGER*4
#  ifdef _I8_
#    define MOLCAS_TO_BLAS_INT
#  endif
#else
#  define BLASINT INTEGER*8
#endif


subroutine drot_(n_,dx,incx_,dy,incy_,c,s)
  implicit none
  integer n_, incx_, incy_
  real*8 dx(*),dy(*),c,s
#ifdef MOLCAS_TO_BLAS_INT
  BLASINT n,  incx,  incy
  n=n_
  incx=incx_
  incy=incy_
  call drot(n,dx,incx,dy,incy,c,s)
#else
  call drot(n_,dx,incx_,dy,incy_,c,s)
#endif
end subroutine

subroutine dswap_(n_,dx,incx_,dy,incy_)
  implicit none
  integer n_, incx_, incy_
  real*8 dx(*), dy(*)
#ifdef MOLCAS_TO_BLAS_INT
  BLASINT n,  incx,  incy
  n=n_
  incx=incx_
  incy=incy_
  call dswap(n,dx,incx,dy,incy)
#else
  call dswap(n_,dx,incx_,dy,incy_)
#endif
end subroutine

subroutine dscal_(n_,da,dx,incx_)
  implicit none
  integer n_, incx_
  real*8 da, dx(*)
#ifdef MOLCAS_TO_BLAS_INT
  BLASINT n,  incx
  n=n_
  incx=incx_
  call dscal(n,da,dx,incx)
#else
  call dscal(n_,da,dx,incx_)
#endif
end subroutine

subroutine dcopy_(n_,dx,incx_,dy,incy_)
  implicit none
  integer n_, incx_, incy_
  real*8 dx(*), dy(*)
#ifdef MOLCAS_TO_BLAS_INT
  BLASINT n,  incx,  incy
  n=n_
  incx=incx_
  incy=incy_
  call dcopy(n,dx,incx,dy,incy)
#else
  call dcopy(n_,dx,incx_,dy,incy_)
#endif
end subroutine

subroutine zcopy_(n_,dx,incx_,dy,incy_)
  implicit none
  integer n_, incx_, incy_
  complex*16 dx(*), dy(*)
#ifdef MOLCAS_TO_BLAS_INT
  BLASINT n,  incx,  incy
  n=n_
  incx=incx_
  incy=incy_
  call zcopy(n,dx,incx,dy,incy)
#else
  call zcopy(n_,dx,incx_,dy,incy_)
#endif
end subroutine

subroutine scopy_(n_,sx,incx_,sy,incy_)
  implicit none
  integer n_, incx_, incy_
  real*4 sx(*), sy(*)
#ifdef MOLCAS_TO_BLAS_INT
  BLASINT n,  incx,  incy
  n=n_
  incx=incx_
  incy=incy_
  call scopy(n,sx,incx,sy,incy)
#else
  call scopy(n_,sx,incx_,sy,incy_)
#endif
end subroutine

subroutine daxpy_(n_,da,dx,incx_,dy,incy_)
  implicit none
  integer n_, incx_, incy_
  real*8 da, dx(*), dy(*)
#ifdef MOLCAS_TO_BLAS_INT
  BLASINT n,  incx,  incy
  n=n_
  incx=incx_
  incy=incy_
  call daxpy(n,da,dx,incx,dy,incy)
#else
  call daxpy(n_,da,dx,incx_,dy,incy_)
#endif
end subroutine

real*8 function ddot_(n_,dx,incx_,dy,incy_)
  implicit none
  integer n_, incx_, incy_
  real*8 dx(*), dy(*)
  real*8 ddot
#ifdef MOLCAS_TO_BLAS_INT
  BLASINT n,  incx,  incy
  n=n_
  incx=incx_
  incy=incy_
  ddot_ = ddot(n,dx,incx,dy,incy)
#else
  ddot_ = ddot(n_,dx,incx_,dy,incy_)
#endif
end function

real*8 function dnrm2_(n_,x,incx_)
  implicit none
  integer n_, incx_
  real*8 x(*)
  real*8, external :: dnrm2
#ifdef MOLCAS_TO_BLAS_INT
  BLASINT n,  incx
  n=n_
  incx=incx_
  dnrm2_ = dnrm2(n,x,incx)
#else
  dnrm2_ = dnrm2(n_,x,incx_)
#endif
end function


real*8 function dznrm2_(n_,x,incx_)
  implicit none
  integer n_, incx_
  complex*16 x(*)
  real*8, external :: dznrm2
#ifdef MOLCAS_TO_BLAS_INT
  BLASINT n,  incx
  n=n_
  incx=incx_
  dznrm2_ = dznrm2(n,x,incx)
#else
  dznrm2_ = dznrm2(n_,x,incx_)
#endif
end function



real*8 function dasum_(n_,dx,incx_)
  implicit none
  integer n_, incx_
  real*8 dx(*)
  real*8 dasum
#ifdef MOLCAS_TO_BLAS_INT
  BLASINT n,  incx
  n=n_
  incx=incx_
  dasum_ = dasum(n,dx,incx)
#else
  dasum_ = dasum(n_,dx,incx_)
#endif
end function

integer function idamax_(n_,dx,incx_)
  implicit none
  integer n_, incx_
  real*8 dx(*)
#ifdef MOLCAS_TO_BLAS_INT
  BLASINT n,  incx
  BLASINT, external :: idamax
  n=n_
  incx=incx_
  idamax_ = idamax(n,dx,incx)
#else
  integer, external :: idamax
  idamax_ = idamax(n_,dx,incx_)
#endif
end function

subroutine dgemm_(transa,transb,m_,n_,k_, &
      &           alpha,a,lda_,b,ldb_,beta,c,ldc_)
  implicit none
  character        transa, transb
  real*8           alpha, beta
  integer          k_, lda_, ldb_, ldc_, m_, n_
  real*8           a(lda_,*), b(ldb_,*), c(ldc_,*)
#ifdef MOLCAS_TO_BLAS_INT
  BLASINT          k,  lda,  ldb,  ldc,  m,  n
#endif
#ifdef _CUDA_BLAS_
  parameter (ncuda=128*128)
  integer*4        k4,lda4,ldb4,ldc4,m4,n4
#endif

  if(m_.eq.0.and.n_.eq.0) return

#ifdef _CUDA_BLAS_
  if(n*m.gt.ncuda) then
    k4=k_
    lda4=lda_
    ldb4=ldb_
    ldc4=ldc_
    m4=m_
    n4=n_
    call cublas_dgemm(transa,transb,
    &     m4,n4,k4,alpha,a,lda4,b,ldb4,beta,c,ldc4)
  else
#endif

#ifdef MOLCAS_TO_BLAS_INT
    m=m_
    n=n_
    k=k_
    lda=lda_
    ldb=ldb_
    ldc=ldc_
    call dgemm(transa,transb, &
        &      m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
#else
    call dgemm(transa,transb,m_,n_,k_, &
        &      alpha,a,lda_,b,ldb_,beta,c,ldc_)
#endif

#ifdef _CUDA_BLAS_
  endif
#endif
end subroutine


subroutine zgemm_(transa,transb,m_,n_,k_, &
      &           alpha,a,lda_,b,ldb_,beta,c,ldc_)
  implicit none
  character        transa, transb
  complex*16       alpha, beta
  integer          k_, lda_, ldb_, ldc_, m_, n_
  complex*16       a(lda_,*), b(ldb_,*), c(ldc_,*)
#ifdef MOLCAS_TO_BLAS_INT
  BLASINT          k,  lda,  ldb,  ldc,  m,  n
#endif
#ifdef _CUDA_BLAS_
  parameter (ncuda=128*128)
  integer*4        k4,lda4,ldb4,ldc4,m4,n4
#endif

  if(m_.eq.0.and.n_.eq.0) return

#ifdef _CUDA_BLAS_
  if(n*m.gt.ncuda) then
    k4=k_
    lda4=lda_
    ldb4=ldb_
    ldc4=ldc_
    m4=m_
    n4=n_
    call cublas_zgemm(transa,transb,
    &     m4,n4,k4,alpha,a,lda4,b,ldb4,beta,c,ldc4)
  else
#endif

#ifdef MOLCAS_TO_BLAS_INT
    m=m_
    n=n_
    k=k_
    lda=lda_
    ldb=ldb_
    ldc=ldc_
    call zgemm(transa,transb, &
        &      m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
#else
    call zgemm(transa,transb,m_,n_,k_, &
        &      alpha,a,lda_,b,ldb_,beta,c,ldc_)
#endif

#ifdef _CUDA_BLAS_
  endif
#endif
end subroutine


subroutine dspmv_ ( uplo,n_,alpha,ap,x,incx_,beta,y,incy_)
  implicit none
  character          uplo
  integer            n_, incx_, incy_
  real*8             alpha, beta
  real*8             ap( * ), x( * ), y( * )
#ifdef MOLCAS_TO_BLAS_INT
  BLASINT            n,  incx,  incy
#endif
#ifdef _CUDA_BLAS_
  parameter (ncuda=128)
  integer*4          n4, incx4, incy4
#endif

#ifdef _CUDA_BLAS_
  if(n.gt.ncuda) then
    n4=n_
    incx4=incx_
    incy4=incy_
    call cublas_dspmv (uplo,n4,alpha,ap,x,incx4,beta,y,incy4)
  else
#endif

#ifdef MOLCAS_TO_BLAS_INT
    n=n_
    incx=incx_
    incy=incy_
    call dspmv ( uplo, n, alpha, ap, x, incx, beta, y, incy )
#else
    call dspmv ( uplo,n_,alpha,ap,x,incx_,beta,y,incy_)
#endif

#ifdef _CUDA_BLAS_
  endif
#endif
end subroutine

subroutine dgemv_(trans,m_,n_,alpha,a,lda_,x,incx_,beta,y,incy_)
  implicit none
  character        trans
  real*8           alpha,beta
  integer          incx_,incy_,lda_,m_,n_
  real*8           a(lda_,*),x(*),y(*)
#ifdef MOLCAS_TO_BLAS_INT
  BLASINT          incx, incy, lda, m, n
#endif
#ifdef _CUDA_BLAS_
  parameter (ncuda=128*128)
  integer*4        lda4,m4,n4,incx4,incy4
#endif

#ifdef _CUDA_BLAS_
  if(n*m.gt.ncuda) then
    m4=m_
    n4=n_
    lda4=lda_
    incx4=incx_
    incy4=incy_
    call cublas_dgemv(trans,m4,n4,alpha,a,lda4,x,incx4,beta,y,incy4)
  else
#endif

#ifdef MOLCAS_TO_BLAS_INT
    m=m_
    n=n_
    lda=lda_
    incx=incx_
    incy=incy_
    call dgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
#else
    call dgemv(trans,m_,n_,alpha,a,lda_,x,incx_,beta,y,incy_)
#endif

#ifdef _CUDA_BLAS_
  endif
#endif
end subroutine
