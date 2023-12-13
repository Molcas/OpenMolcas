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
! dasum_
! daxpy_
! dcopy_
! ddot_
! dgemm_
! dgemv_
! dnrm2_
! drot_
! dscal_
! dspmv_
! dswap_
! dznrm2_
! idamax_
! scopy_
! zaxpy_
! zcopy_
! zdscal_
! zgemm_
! zscal_

! Specify if integer conversion will be needed.
! (real conversion is not implemented yet)
#if defined (LINALG_I4) && defined (_I8_)
# define MOLCAS_TO_BLAS_INT
# define _BLAS_INT_use_ use Definitions, only: BLASInt
#else
# define _BLAS_INT_use_
#endif
#ifdef _CUDA_BLAS_
# define _CUDA_INT_use_ use Definitions, only: CUDAInt
#else
# define _CUDA_INT_use_
#endif

#include "intent.fh"

function dasum_(n_,dx,incx_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  real(kind=BLASR8) :: dasum_
  integer(kind=iwp), intent(in) :: n_, incx_
  real(kind=BLASR8), intent(in) :: dx(*)
  real(kind=BLASR8), external :: dasum
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: n, incx
  n = int(n_,kind=BLASInt)
  incx = int(incx_,kind=BLASInt)
  dasum_ = dasum(n,dx,incx)
# else
  dasum_ = dasum(n_,dx,incx_)
# endif
end function dasum_

subroutine daxpy_(n_,da,dx,incx_,dy,incy_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  integer(kind=iwp), intent(in) :: n_, incx_, incy_
  real(kind=BLASR8), intent(in) :: da, dx(*)
  real(kind=BLASR8), intent(inout) :: dy(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: n, incx, incy
  n = int(n_,kind=BLASInt)
  incx = int(incx_,kind=BLASInt)
  incy = int(incy_,kind=BLASInt)
  call daxpy(n,da,dx,incx,dy,incy)
# else
  call daxpy(n_,da,dx,incx_,dy,incy_)
# endif
end subroutine daxpy_

subroutine dcopy_(n_,dx,incx_,dy,incy_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  integer(kind=iwp), intent(in) :: n_, incx_, incy_
  real(kind=BLASR8), intent(in) :: dx(*)
  real(kind=BLASR8), intent(_OUT_) :: dy(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: n, incx, incy
  n = int(n_,kind=BLASInt)
  incx = int(incx_,kind=BLASInt)
  incy = int(incy_,kind=BLASInt)
  call dcopy(n,dx,incx,dy,incy)
# else
  call dcopy(n_,dx,incx_,dy,incy_)
# endif
end subroutine dcopy_

function ddot_(n_,dx,incx_,dy,incy_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  real(kind=BLASR8) :: ddot_
  integer(kind=iwp), intent(in) :: n_, incx_, incy_
  real(kind=BLASR8), intent(in) :: dx(*), dy(*)
  real(kind=BLASR8), external :: ddot
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: n, incx, incy
  n = int(n_,kind=BLASInt)
  incx = int(incx_,kind=BLASInt)
  incy = int(incy_,kind=BLASInt)
  ddot_ = ddot(n,dx,incx,dy,incy)
# else
  ddot_ = ddot(n_,dx,incx_,dy,incy_)
# endif
end function ddot_

subroutine dgemm_(transa,transb,m_,n_,k_,alpha,a,lda_,b,ldb_,beta,c,ldc_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  _CUDA_INT_use_
  implicit none
  character, intent(in) :: transa, transb
  integer(kind=iwp), intent(in) :: m_, n_, k_, lda_, ldb_, ldc_
  real(kind=BLASR8), intent(in) :: alpha, a(lda_,*), b(ldb_,*), beta
  real(kind=BLASR8), intent(inout) :: c(ldc_,*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: k, lda, ldb, ldc, m, n
# endif
# ifdef _CUDA_BLAS_
  integer(kind=CUDAInt) :: k4, lda4, ldb4, ldc4, m4, n4
  integer(kind=iwp), parameter :: ncuda = 128*128
# endif

  if ((m_ == 0) .and. (n_ == 0)) return

# ifdef _CUDA_BLAS_
  if (n_*m_ > ncuda) then
    m4 = int(m_,kind=CUDAInt)
    n4 = int(n_,kind=CUDAInt)
    k4 = int(k_,kind=CUDAInt)
    lda4 = int(lda_,kind=CUDAInt)
    ldb4 = int(ldb_,kind=CUDAInt)
    ldc4 = int(ldc_,kind=CUDAInt)
    call cublas_dgemm(transa,transb,m4,n4,k4,alpha,a,lda4,b,ldb4,beta,c,ldc4)
  else
# endif

#   ifdef MOLCAS_TO_BLAS_INT
    m = int(m_,kind=BLASInt)
    n = int(n_,kind=BLASInt)
    k = int(k_,kind=BLASInt)
    lda = int(lda_,kind=BLASInt)
    ldb = int(ldb_,kind=BLASInt)
    ldc = int(ldc_,kind=BLASInt)
    call dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
#   else
    call dgemm(transa,transb,m_,n_,k_,alpha,a,lda_,b,ldb_,beta,c,ldc_)
#   endif

# ifdef _CUDA_BLAS_
  end if
# endif
end subroutine dgemm_

subroutine dgemv_(trans,m_,n_,alpha,a,lda_,x,incx_,beta,y,incy_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  _CUDA_INT_use_
  implicit none
  character, intent(in) :: trans
  integer(kind=iwp), intent(in) :: m_, n_, lda_, incx_, incy_
  real(kind=BLASR8), intent(in) :: alpha, a(lda_,*), x(*), beta
  real(kind=BLASR8), intent(inout) :: y(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: incx, incy, lda, m, n
# endif
# ifdef _CUDA_BLAS_
  integer(kind=CUDAInt) :: lda4, m4, n4, incx4, incy4
  integer(kind=iwp), parameter :: ncuda = 128*128
# endif

# ifdef _CUDA_BLAS_
  if (n_*m_ > ncuda) then
    m4 = int(m_,kind=CUDAInt)
    n4 = int(n_,kind=CUDAInt)
    lda4 = int(lda_,kind=CUDAInt)
    incx4 = int(incx_,kind=CUDAInt)
    incy4 = int(incy_,kind=CUDAInt)
    call cublas_dgemv(trans,m4,n4,alpha,a,lda4,x,incx4,beta,y,incy4)
  else
# endif

#   ifdef MOLCAS_TO_BLAS_INT
    m = int(m_,kind=BLASInt)
    n = int(n_,kind=BLASInt)
    lda = int(lda_,kind=BLASInt)
    incx = int(incx_,kind=BLASInt)
    incy = int(incy_,kind=BLASInt)
    call dgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
#   else
    call dgemv(trans,m_,n_,alpha,a,lda_,x,incx_,beta,y,incy_)
#   endif

# ifdef _CUDA_BLAS_
  end if
# endif
end subroutine dgemv_

function dnrm2_(n_,x,incx_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  real(kind=BLASR8) :: dnrm2_
  integer(kind=iwp), intent(in) :: n_, incx_
  real(kind=BLASR8), intent(in) :: x(*)
  real(kind=BLASR8), external :: dnrm2
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: n, incx
  n = int(n_,kind=BLASInt)
  incx = int(incx_,kind=BLASInt)
  dnrm2_ = dnrm2(n,x,incx)
# else
  dnrm2_ = dnrm2(n_,x,incx_)
# endif
end function dnrm2_

subroutine drot_(n_,dx,incx_,dy,incy_,c,s)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  integer(kind=iwp), intent(in) :: n_, incx_, incy_
  real(kind=BLASR8), intent(inout) :: dx(*), dy(*)
  real(kind=BLASR8), intent(in) :: c, s
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: n, incx, incy
  n = int(n_,kind=BLASInt)
  incx = int(incx_,kind=BLASInt)
  incy = int(incy_,kind=BLASInt)
  call drot(n,dx,incx,dy,incy,c,s)
# else
  call drot(n_,dx,incx_,dy,incy_,c,s)
# endif
end subroutine drot_

subroutine dscal_(n_,da,dx,incx_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  integer(kind=iwp), intent(in) :: n_, incx_
  real(kind=BLASR8), intent(in) :: da
  real(kind=BLASR8), intent(inout) :: dx(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: n, incx
  n = int(n_,kind=BLASInt)
  incx = int(incx_,kind=BLASInt)
  call dscal(n,da,dx,incx)
# else
  call dscal(n_,da,dx,incx_)
# endif
end subroutine dscal_

subroutine dspmv_(uplo,n_,alpha,ap,x,incx_,beta,y,incy_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  _CUDA_INT_use_
  implicit none
  character, intent(in) :: uplo
  integer(kind=iwp), intent(in) :: n_, incx_, incy_
  real(kind=BLASR8), intent(in) :: alpha, ap(*), x(*), beta
  real(kind=BLASR8), intent(inout) :: y(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: n, incx, incy
# endif
# ifdef _CUDA_BLAS_
  integer(kind=CUDAInt) :: n4, incx4, incy4
  integer(kind=iwp), parameter :: ncuda = 128
# endif

# ifdef _CUDA_BLAS_
  if (n_ > ncuda) then
    n4 = int(n_,kind=CUDAInt)
    incx4 = int(incx_,kind=CUDAInt)
    incy4 = int(incy_,kind=CUDAInt)
    call cublas_dspmv(uplo,n4,alpha,ap,x,incx4,beta,y,incy4)
  else
# endif

#   ifdef MOLCAS_TO_BLAS_INT
    n = int(n_,kind=BLASInt)
    incx = int(incx_,kind=BLASInt)
    incy = int(incy_,kind=BLASInt)
    call dspmv(uplo,n,alpha,ap,x,incx,beta,y,incy)
#   else
    call dspmv(uplo,n_,alpha,ap,x,incx_,beta,y,incy_)
#   endif

# ifdef _CUDA_BLAS_
  end if
# endif
end subroutine dspmv_

subroutine dswap_(n_,dx,incx_,dy,incy_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  integer(kind=iwp), intent(in) :: n_, incx_, incy_
  real(kind=BLASR8), intent(inout) :: dx(*), dy(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: n, incx, incy
  n = int(n_,kind=BLASInt)
  incx = int(incx_,kind=BLASInt)
  incy = int(incy_,kind=BLASInt)
  call dswap(n,dx,incx,dy,incy)
# else
  call dswap(n_,dx,incx_,dy,incy_)
# endif
end subroutine dswap_

function dznrm2_(n_,x,incx_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  real(kind=BLASR8) :: dznrm2_
  integer(kind=iwp), intent(in) :: n_, incx_
  complex(kind=BLASR8), intent(in) :: x(*)
  real(kind=BLASR8), external :: dznrm2
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: n, incx
  n = int(n_,kind=BLASInt)
  incx = int(incx_,kind=BLASInt)
  dznrm2_ = dznrm2(n,x,incx)
# else
  dznrm2_ = dznrm2(n_,x,incx_)
# endif
end function dznrm2_

function idamax_(n_,dx,incx_)
  use Definitions, only: BLASR8, iwp, BLASInt
  implicit none
  integer(kind=iwp) :: idamax_
  integer(kind=iwp), intent(in) :: n_, incx_
  real(kind=BLASR8), intent(in) :: dx(*)
  integer(kind=BLASInt), external :: idamax
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: n, incx
  n = int(n_,kind=BLASInt)
  incx = int(incx_,kind=BLASInt)
  idamax_ = idamax(n,dx,incx)
# else
  idamax_ = idamax(n_,dx,incx_)
# endif
end function idamax_

subroutine scopy_(n_,sx,incx_,sy,incy_)
  use Definitions, only: BLASR4, iwp
  _BLAS_INT_use_
  implicit none
  integer(kind=iwp), intent(in) :: n_, incx_, incy_
  real(kind=BLASR4), intent(in) :: sx(*)
  real(kind=BLASR4), intent(_OUT_) :: sy(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: n, incx, incy
  n = int(n_,kind=BLASInt)
  incx = int(incx_,kind=BLASInt)
  incy = int(incy_,kind=BLASInt)
  call scopy(n,sx,incx,sy,incy)
# else
  call scopy(n_,sx,incx_,sy,incy_)
# endif
end subroutine scopy_

subroutine zaxpy_(n_,da,dx,incx_,dy,incy_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  integer(kind=iwp), intent(in) :: n_, incx_, incy_
  complex(kind=BLASR8), intent(in) :: da, dx(*)
  complex(kind=BLASR8), intent(out) :: dy(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: n, incx, incy
  n = int(n_,kind=BLASInt)
  incx = int(incx_,kind=BLASInt)
  incy = int(incy_,kind=BLASInt)
  call zaxpy(n,da,dx,incx,dy,incy)
# else
  call zaxpy(n_,da,dx,incx_,dy,incy_)
# endif
end subroutine zaxpy_

subroutine zcopy_(n_,dx,incx_,dy,incy_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  integer(kind=iwp), intent(in) :: n_, incx_, incy_
  complex(kind=BLASR8), intent(in) :: dx(*)
  complex(kind=BLASR8), intent(_OUT_) :: dy(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: n, incx, incy
  n = int(n_,kind=BLASInt)
  incx = int(incx_,kind=BLASInt)
  incy = int(incy_,kind=BLASInt)
  call zcopy(n,dx,incx,dy,incy)
# else
  call zcopy(n_,dx,incx_,dy,incy_)
# endif
end subroutine zcopy_

subroutine zdscal_(n_,da,dx,incx_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  integer(kind=iwp), intent(in) :: n_, incx_
  real(kind=BLASR8), intent(in) :: da
  complex(kind=BLASR8), intent(inout) :: dx(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: n, incx
  n = int(n_,kind=BLASInt)
  incx = int(incx_,kind=BLASInt)
  call zdscal(n,da,dx,incx)
# else
  call zdscal(n_,da,dx,incx_)
# endif
end subroutine zdscal_

subroutine zgemm_(transa,transb,m_,n_,k_,alpha,a,lda_,b,ldb_,beta,c,ldc_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  _CUDA_INT_use_
  implicit none
  character, intent(in) :: transa, transb
  integer(kind=iwp), intent(in) :: m_, n_, k_, lda_, ldb_, ldc_
  complex(kind=BLASR8), intent(in) :: alpha, a(lda_,*), b(ldb_,*), beta
  complex(kind=BLASR8), intent(inout) :: c(ldc_,*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: k, lda, ldb, ldc, m, n
# endif
# ifdef _CUDA_BLAS_
  integer(kind=CUDAInt) :: k4, lda4, ldb4, ldc4, m4, n4
  integer(kind=iwp), parameter :: ncuda = 128*128
# endif

  if ((m_ == 0) .and. (n_ == 0)) return

# ifdef _CUDA_BLAS_
  if (n_*m_ > ncuda) then
    m4 = int(m_,kind=CUDAInt)
    n4 = int(n_,kind=CUDAInt)
    k4 = int(k_,kind=CUDAInt)
    lda4 = int(lda_,kind=CUDAInt)
    ldb4 = int(ldb_,kind=CUDAInt)
    ldc4 = int(ldc_,kind=CUDAInt)
    call cublas_zgemm(transa,transb,m4,n4,k4,alpha,a,lda4,b,ldb4,beta,c,ldc4)
  else
# endif

#   ifdef MOLCAS_TO_BLAS_INT
    m = int(m_,kind=BLASInt)
    n = int(n_,kind=BLASInt)
    k = int(k_,kind=BLASInt)
    lda = int(lda_,kind=BLASInt)
    ldb = int(ldb_,kind=BLASInt)
    ldc = int(ldc_,kind=BLASInt)
    call zgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
#   else
    call zgemm(transa,transb,m_,n_,k_,alpha,a,lda_,b,ldb_,beta,c,ldc_)
#   endif

# ifdef _CUDA_BLAS_
  end if
# endif
end subroutine zgemm_

subroutine zscal_(n_,da,dx,incx_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  integer(kind=iwp), intent(in) :: n_, incx_
  complex(kind=BLASR8), intent(in) :: da
  complex(kind=BLASR8), intent(inout) :: dx(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: n, incx
  n = int(n_,kind=BLASInt)
  incx = int(incx_,kind=BLASInt)
  call zscal(n,da,dx,incx)
# else
  call zscal(n_,da,dx,incx_)
# endif
end subroutine zscal_
