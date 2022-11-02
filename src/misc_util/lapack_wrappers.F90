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
! Molcas LAPACK wrappers
! These wrappers handle lp64/ilp64 interfaces and various other
! requirements related to matrix size on GPUs. Any call to LAPACK in
! Molcas should use these wrappers.

! List of LAPACK wrappers currently implemented:
! dopmtr_
! dspgv_
! dsptrd_
! dstevr_
! dgetrs_
! dspev_
! dgees_
! dgeev_
! dgesvd_
! dgetrf_
! dgesv_
! dsyev_
! dsyevr_
! dgetri_
! zhpev_
! dsygv_
! dpotrf_
! dgels_
! dsytrd_
! dposv_
! dlamch_
! ilaenv_
! dsterf_
! dsteqr_
! dlascl_
! dorgtr_
! zgesvd_

! Specify if integer conversion will be needed.
#if defined(LINALG_I4) && defined(_I8_)
# define MOLCAS_TO_BLAS_INT
# define _BLAS_INT_use_ use Definitions, only: BLASInt
# define _BLAS_INT_stdalloc_ \
  use stdalloc, only: mma_allocate, mma_deallocate
#else
# define _BLAS_INT_use_
# define _BLAS_INT_stdalloc_ !
#endif

! For procedures known to raise floating point exceptions in the test suite,
! disable exception trapping locally: three pieces of code are needed
#ifdef _FPE_TRAP_
  ! can't use "only" in IEEE_exceptions, because of line length
# define _FPE_TRAP_use_ \
  use, intrinsic :: IEEE_exceptions; \
  use Definitions, only: DI => DefInt
# define _FPE_TRAP_init_ \
  type(IEEE_Status_Type) :: IEEE_Status; \
  call IEEE_Get_Status(IEEE_Status); \
  call IEEE_Set_Halting_Mode(IEEE_Usual,.false._DI)
# define _FPE_TRAP_end_ \
  call IEEE_Set_Status(IEEE_Status)
#else
# define _FPE_TRAP_use_ !
# define _FPE_TRAP_init_ !
# define _FPE_TRAP_end_ !
#endif

#include "macros.fh"

subroutine dopmtr_(side,uplo,trans,m_,n_,ap,tau,c,ldc_,work,info_)
  use Definitions, only: wp, iwp
  _BLAS_INT_use_
  implicit none
  character :: side, uplo, trans
  integer(kind=iwp) :: m_, n_, ldc_, info_
  real(kind=wp) :: ap(*), tau(*), c(ldc_,*), work(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: info, ldc, m, n
  m = int(m_,kind=BLASInt)
  n = int(n_,kind=BLASInt)
  ldc = int(ldc_,kind=BLASInt)
  call dopmtr(side,uplo,trans,m,n,ap,tau,c,ldc,work,info)
  info_ = info
# else
  call dopmtr(side,uplo,trans,m_,n_,ap,tau,c,ldc_,work,info_)
# endif
end subroutine dopmtr_

subroutine dspgv_(itype_,jobz,uplo,n_,ap,bp,w,z,ldz_,work,info_)
  use Definitions, only: wp, iwp
  _BLAS_INT_use_
  implicit none
  integer(kind=iwp) :: itype_, n_, ldz_, info_
  character :: jobz, uplo
  real(kind=wp) :: ap(*), bp(*), w(*), z(ldz_,*), work(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: info, itype, ldz, n
  itype = int(itype_,kind=BLASInt)
  n = int(n_,kind=BLASInt)
  ldz = int(ldz_,kind=BLASInt)
  call dspgv(itype,jobz,uplo,n,ap,bp,w,z,ldz,work,info)
  info_ = info
# else
  call dspgv(itype_,jobz,uplo,n_,ap,bp,w,z,ldz_,work,info_)
# endif
end subroutine dspgv_

subroutine dsptrd_(uplo,n_,ap,d,e,tau,info_)
  use Definitions, only: wp, iwp
  _BLAS_INT_use_
  implicit none
  character :: uplo
  integer(kind=iwp) :: n_, info_
  real(kind=wp) :: ap(*), d(*), e(*), tau(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: n, info
  n = int(n_,kind=BLASInt)
  call dsptrd(uplo,n,ap,d,e,tau,info)
  info_ = info
# else
  call dsptrd(uplo,n_,ap,d,e,tau,info_)
# endif
end subroutine dsptrd_

subroutine dstevr_(jobz,rng,n_,d,e,vl,vu,il_,iu_,abstol,m_,w,z,ldz_,isuppz_,work,lwork_,iwork_,liwork_,info_)
  use Definitions, only: wp, iwp
  _BLAS_INT_use_
  _BLAS_INT_stdalloc_
  _FPE_TRAP_use_
  implicit none
  character :: jobz, rng
  integer(kind=iwp) :: n_, il_, iu_, m_, ldz_, isuppz_(*), lwork_, iwork_(*), liwork_, info_
  real(kind=wp) :: d(*), e(*), vl, vu, abstol, w(*), z(ldz_,*), work(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=iwp) :: i
  integer(kind=BLASInt) :: il, info, iu, ldz, liwork, lwork, m, n
  integer(kind=BLASInt), allocatable :: isuppz(:), iwork(:)
  _FPE_TRAP_init_
  unused_var(lwork_)
  unused_var(iwork_(1))
  n = int(n_,kind=BLASInt)
  il = int(il_,kind=BLASInt)
  iu = int(iu_,kind=BLASInt)
  m = int(m_,kind=BLASInt)
  ldz = int(ldz_,kind=BLASInt)
  liwork = int(liwork_,kind=BLASInt)
  i=2*max(1,m_)
  call mma_allocate(isuppz,i,label='isuppz')
  call mma_allocate(iwork,liwork_,label='iwork')
  call dstevr(jobz,rng,n,d,e,vl,vu,il,iu,abstol,m,w,z,ldz,isuppz,work,lwork,iwork,liwork,info)
  call mma_deallocate(iwork)
  isuppz_(1:i) = isuppz(1:i)
  call mma_deallocate(isuppz)
  info_ = info
# else
  _FPE_TRAP_init_
  call dstevr(jobz,rng,n_,d,e,vl,vu,il_,iu_,abstol,m_,w,z,ldz_,isuppz_,work,lwork_,iwork_,liwork_,info_)
# endif
  _FPE_TRAP_end_
end subroutine dstevr_

subroutine dgetrs_(trans,n_,nrhs_,a,lda_,ipiv_,b,ldb_,info_)
  use Definitions, only: wp, iwp
  _BLAS_INT_use_
  _BLAS_INT_stdalloc_
  implicit none
  character :: trans
  integer(kind=iwp) :: n_, nrhs_, lda_, ipiv_(*), ldb_, info_
  real(kind=wp) :: a(lda_,*), b(ldb_,*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: i, info, lda, ldb, n, nrhs
  integer(kind=BLASInt), allocatable :: ipiv(:)
  n = int(n_,kind=BLASInt)
  nrhs = int(nrhs_,kind=BLASInt)
  lda = int(lda_,kind=BLASInt)
  call mma_allocate(ipiv,n_,label='ipiv')
  do i=1,n
    ipiv(i) = int(ipiv_(i),kind=BLASInt)
  end do
  ldb = int(ldb_,kind=BLASInt)
  call dgetrs(trans,n,nrhs,a,lda,ipiv,b,ldb,info)
  call mma_deallocate(ipiv)
  info_ = info
# else
  call dgetrs(trans,n_,nrhs_,a,lda_,ipiv_,b,ldb_,info_)
# endif
end subroutine dgetrs_

subroutine dspev_(jobz,uplo,n_,ap,w,z,ldz_,work,info_)
  use Definitions, only: wp, iwp
  _BLAS_INT_use_
  implicit none
  character :: jobz, uplo
  integer(kind=iwp) :: n_, ldz_, info_
  real(kind=wp) :: ap(*), w(*), z(ldz_,*), work(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: info, ldz, n
  n = int(n_,kind=BLASInt)
  ldz = int(ldz_,kind=BLASInt)
  call dspev(jobz,uplo,n,ap,w,z,ldz,work,info)
  info_ = info
# else
  call dspev(jobz,uplo,n_,ap,w,z,ldz_,work,info_)
# endif
end subroutine dspev_

subroutine dgees_(jobvs,sort,slct,n_,a,lda_,sdim_,wr,wi,vs,ldvs_,work,lwork_,bwork,info_)
  use Definitions, only: wp, iwp
  _BLAS_INT_use_
  implicit none
  character :: jobvs, sort
  integer(kind=iwp) :: n_, lda_, sdim_, ldvs_, lwork_, info_
  real(kind=wp) :: a(lda_,*), wr(*), wi(*), vs(ldvs_,*), work(*)
  logical(kind=iwp), external :: slct
  logical(kind=iwp) :: bwork(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: n, lda, sdim, ldvs, lwork, info
  n = int(n_,kind=BLASInt)
  lda = int(lda_,kind=BLASInt)
  ldvs = int(ldvs_,kind=BLASInt)
  lwork = int(lwork_,kind=BLASInt)
  call dgees(jobvs,sort,slct,n,a,lda,sdim,wr,wi,vs,ldvs,work,lwork,bwork,info)
  info_ = info
  sdim_ = sdim
# else
  call dgees(jobvs,sort,slct,n_,a,lda_,sdim_,wr,wi,vs,ldvs_,work,lwork_,bwork,info_)
# endif
end subroutine dgees_

subroutine dgeev_(jobvl,jobvr,n_,a,lda_,wr,wi,vl,ldvl_,vr,ldvr_,work,lwork_,info_)
  use Definitions, only: wp, iwp
  _BLAS_INT_use_
  implicit none
  character :: jobvl, jobvr
  integer(kind=iwp) :: n_, lda_, ldvl_, ldvr_, lwork_, info_
  real(kind=wp) :: a(lda_,*), wr(*), wi(*), vl(ldvl_,*), vr(ldvr_,*), work(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: info, lda, ldvl, ldvr, lwork, n
  n = int(n_,kind=BLASInt)
  lda = int(lda_,kind=BLASInt)
  ldvl = int(ldvl_,kind=BLASInt)
  ldvr = int(ldvr_,kind=BLASInt)
  lwork = int(lwork_,kind=BLASInt)
  call dgeev(jobvl,jobvr,n,a,lda,wr,wi,vl,ldvl,vr,ldvr,work,lwork,info)
  info_ = info
# else
  call dgeev(jobvl,jobvr,n_,a,lda_,wr,wi,vl,ldvl_,vr,ldvr_,work,lwork_,info_)
# endif
end subroutine dgeev_

subroutine dgesvd_(jobu,jobvt,m_,n_,a,lda_,s,u,ldu_,vt,ldvt_,work,lwork_,info_)
  use Definitions, only: wp, iwp
  _BLAS_INT_use_
  implicit none
  character :: jobu, jobvt
  integer(kind=iwp) :: m_, n_, lda_, ldu_, ldvt_, lwork_, info_
  real(kind=wp) :: a(lda_,*), s(*), u(ldu_,*), vt(ldvt_,*), work(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: info, lda, ldu, ldvt, lwork, m, n
  m = int(m_,kind=BLASInt)
  n = int(n_,kind=BLASInt)
  lda = int(lda_,kind=BLASInt)
  ldu = int(ldu_,kind=BLASInt)
  ldvt = int(ldvt_,kind=BLASInt)
  lwork = int(lwork_,kind=BLASInt)
  call dgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,info)
  info_ = info
# else
  call dgesvd(jobu,jobvt,m_,n_,a,lda_,s,u,ldu_,vt,ldvt_,work,lwork_,info_)
# endif
end subroutine dgesvd_

subroutine dgetrf_(m_,n_,a,lda_,ipiv_,info_)
  use Definitions, only: wp, iwp
  _BLAS_INT_use_
  _BLAS_INT_stdalloc_
  implicit none
  integer(kind=iwp) :: m_, n_, lda_, ipiv_(*), info_
  real(kind=wp) :: a(lda_,*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: info, lda, m, n
  integer(kind=BLASInt), allocatable :: ipiv(:)
  m = int(m_,kind=BLASInt)
  n = int(n_,kind=BLASInt)
  lda = int(lda_,kind=BLASInt)
  call mma_allocate(ipiv,n_,label='ipiv')
  call dgetrf(m,n,a,lda,ipiv,info)
  ipiv_(1:n_) = ipiv
  call mma_deallocate(ipiv)
  info_ = info
# else
  call dgetrf(m_,n_,a,lda_,ipiv_,info_)
# endif
end subroutine dgetrf_

subroutine dgesv_(n_,nrhs_,a,lda_,ipiv_,b,ldb_,info_)
  use Definitions, only: wp, iwp
  _BLAS_INT_use_
  _BLAS_INT_stdalloc_
  implicit none
  integer(kind=iwp) :: n_, nrhs_, lda_, ipiv_(*), ldb_, info_
  real(kind=wp) :: a(lda_,*), b(ldb_,*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=iwp) :: i
  integer(kind=BLASInt) :: info, lda, ldb, n, nrhs
  integer(kind=BLASInt), allocatable :: ipiv(:)
  n = int(n_,kind=BLASInt)
  nrhs = int(nrhs_,kind=BLASInt)
  lda = int(lda_,kind=BLASInt)
  call mma_allocate(ipiv,n_,label='ipiv')
  do i=1,n
    ipiv(i) = int(ipiv_(i),kind=BLASInt)
  end do
  ldb = int(ldb_,kind=BLASInt)
  call dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
  call mma_deallocate(ipiv)
  info_ = info
# else
  call dgesv(n_,nrhs_,a,lda_,ipiv_,b,ldb_,info_)
# endif
end subroutine dgesv_

subroutine dsyev_(jobz,uplo,n_,a,lda_,w,work,lwork_,info_)
  use Definitions, only: wp, iwp
  _BLAS_INT_use_
  implicit none
  character :: jobz, uplo
  integer(kind=iwp) :: n_, lda_, lwork_, info_
  real(kind=wp) :: a(lda_,*), w(*), work(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: info, lda, lwork, n
  n = int(n_,kind=BLASInt)
  lda = int(lda_,kind=BLASInt)
  lwork = int(lwork_,kind=BLASInt)
  call dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
  info_ = info
# else
  call dsyev(jobz,uplo,n_,a,lda_,w,work,lwork_,info_)
# endif
end subroutine dsyev_

subroutine dsyevr_(jobz,rng,uplo,n_,a,lda_,vl,vu,il_,iu_,abstol,m_,w,z,ldz_,isuppz_,work,lwork_,iwork_,liwork_,info_)
  use Definitions, only: wp, iwp
  _BLAS_INT_use_
  _BLAS_INT_stdalloc_
  implicit none
  character :: jobz, rng, uplo
  integer(kind=iwp) :: n_, lda_, il_, iu_, m_, ldz_, isuppz_(*), lwork_, iwork_(*), liwork_, info_
  real(kind=wp) :: a(lda_,*), vl, vu, abstol, w(*), z(ldz_,*), work(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=iwp) :: i
  integer(kind=BLASInt) :: il, info, iu, lda, ldz, liwork, lwork, m, n
  integer(kind=BLASInt), allocatable :: isuppz(:), iwork(:)
  n = int(n_,kind=BLASInt)
  lda = int(lda_,kind=BLASInt)
  il = int(il_,kind=BLASInt)
  iu = int(iu_,kind=BLASInt)
  m = int(m_,kind=BLASInt)
  ldz = int(ldz_,kind=BLASInt)
  lwork = int(lwork_,kind=BLASInt)
  liwork = int(liwork_,kind=BLASInt)
  i=2*max(1,m_)
  call mma_allocate(isuppz,i,label='isuppz')
  call mma_allocate(iwork,max(1,liwork_),label='iwork')
  call dsyevr(jobz,rng,uplo,n,a,lda,vl,vu,il,iu,abstol,m,w,z,ldz,isuppz,work,lwork,iwork,liwork,info)
  isuppz_(1:i) = isuppz
  call mma_deallocate(isuppz)
  iwork_(1) = iwork(1)
  call mma_deallocate(iwork)
  info_ = info
# else
  call dsyevr(jobz,rng,uplo,n_,a,lda_,vl,vu,il_,iu_,abstol,m_,w,z,ldz_,isuppz_,work,lwork_,iwork_,liwork_,info_)
# endif
end subroutine dsyevr_

subroutine dgetri_(n_,a,lda_,ipiv_,work,lwork_,info_)
  use Definitions, only: wp, iwp
  _BLAS_INT_use_
  _BLAS_INT_stdalloc_
  implicit none
  integer(kind=iwp) :: n_, lda_, ipiv_(*), lwork_, info_
  real(kind=wp) :: a(lda_,*), work(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=iwp) :: i
  integer(kind=BLASInt) :: info, lda, lwork, n
  integer(kind=BLASInt), allocatable :: ipiv(:)
  n = int(n_,kind=BLASInt)
  lda = int(lda_,kind=BLASInt)
  call mma_allocate(ipiv,n_,label='ipiv')
  do i=1,n
    ipiv(i) = int(ipiv_(i),kind=BLASInt)
  end do
  lwork = int(lwork_,kind=BLASInt)
  call dgetri(n,a,lda,ipiv,work,lwork,info)
  call mma_deallocate(ipiv)
  info_ = info
# else
  call dgetri(n_,a,lda_,ipiv_,work,lwork_,info_)
# endif
end subroutine dgetri_

subroutine zhpev_(jobz,uplo,n_,ap,w,z,ldz_,work,rwork,info_)
  use Definitions, only: wp, iwp
  _BLAS_INT_use_
  implicit none
  character :: jobz, uplo
  integer(kind=iwp) :: n_, ldz_, info_
  complex(kind=wp) :: ap(*), z(ldz_,*), work(*)
  real(kind=wp) :: w(*), rwork(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: info, ldz, n
  n = int(n_,kind=BLASInt)
  ldz = int(ldz_,kind=BLASInt)
  call zhpev(jobz,uplo,n,ap,w,z,ldz,work,rwork,info)
  info_ = info
# else
  call zhpev(jobz,uplo,n_,ap,w,z,ldz_,work,rwork,info_)
# endif
end subroutine zhpev_

subroutine dsygv_(itype_,jobz,uplo,n_,a,lda_,b,ldb_,w,work,lwork_,info_)
  use Definitions, only: wp, iwp
  _BLAS_INT_use_
  implicit none
  integer(kind=iwp) :: itype_, n_, lda_, ldb_, lwork_, info_
  character :: jobz, uplo
  real(kind=wp) :: a(lda_,*), b(ldb_,*), w(*), work(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: info, itype, lda, ldb, lwork, n
  itype = int(itype_,kind=BLASInt)
  n = int(n_,kind=BLASInt)
  lda = int(lda_,kind=BLASInt)
  ldb = int(ldb_,kind=BLASInt)
  lwork = int(lwork_,kind=BLASInt)
  call dsygv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,info)
  info_ = info
# else
  call dsygv(itype_,jobz,uplo,n_,a,lda_,b,ldb_,w,work,lwork_,info_)
# endif
end subroutine dsygv_

subroutine dpotrf_(uplo,n_,a,lda_,info_)
  use Definitions, only: wp, iwp
  _BLAS_INT_use_
  implicit none
  character :: uplo
  integer(kind=iwp) :: n_, lda_, info_
  real(kind=wp) :: a(lda_,*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: info, lda, n
  n = int(n_,kind=BLASInt)
  lda = int(lda_,kind=BLASInt)
  call dpotrf(uplo,n,a,lda,info)
  info_ = info
# else
  call dpotrf(uplo,n_,a,lda_,info_)
# endif
end subroutine dpotrf_

subroutine dgels_(trans,m_,n_,nrhs_,a,lda_,b,ldb_,work,lwork_,info_)
  use Definitions, only: wp, iwp
  _BLAS_INT_use_
  implicit none
  character :: trans
  integer(kind=iwp) :: m_, n_, nrhs_, lda_, ldb_, lwork_, info_
  real(kind=wp) :: a(lda_,*), b(ldb_,*), work(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: info, lda, ldb, lwork, m, n, nrhs
  m = int(m_,kind=BLASInt)
  n = int(n_,kind=BLASInt)
  nrhs = int(nrhs_,kind=BLASInt)
  lda = int(lda_,kind=BLASInt)
  ldb = int(ldb_,kind=BLASInt)
  lwork = int(lwork_,kind=BLASInt)
  call dgels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
  info_ = info
# else
  call dgels(trans,m_,n_,nrhs_,a,lda_,b,ldb_,work,lwork_,info_)
# endif
end subroutine dgels_

subroutine dsytrd_(uplo,n_,a,lda_,d,e,tau,work,lwork_,info_)
  use Definitions, only: wp, iwp
  _BLAS_INT_use_
  implicit none
  character :: uplo
  integer(kind=iwp) :: n_, lda_, lwork_, info_
  real(kind=wp) :: a(lda_,*), d(*), e(*), tau(*), work(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: info, lda, lwork, n
  n = int(n_,kind=BLASInt)
  lda = int(lda_,kind=BLASInt)
  lwork = int(lwork_,kind=BLASInt)
  call dsytrd(uplo,n,a,lda,d,e,tau,work,lwork,info)
  info_ = info
# else
  call dsytrd(uplo,n_,a,lda_,d,e,tau,work,lwork_,info_)
# endif
end subroutine dsytrd_

subroutine dposv_(uplo,n_,nrhs_,a,lda_,b,ldb_,info_)
  use Definitions, only: wp, iwp
  _BLAS_INT_use_
  implicit none
  character :: uplo
  integer(kind=iwp) :: n_, nrhs_, lda_, ldb_, info_
  real(kind=wp) :: a(lda_,*), b(ldb_,*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: info, lda, ldb, n, nrhs
  n = int(n_,kind=BLASInt)
  nrhs = int(nrhs_,kind=BLASInt)
  lda = int(lda_,kind=BLASInt)
  ldb = int(ldb_,kind=BLASInt)
  call dposv(uplo,n,nrhs,a,lda,b,ldb,info)
  info_ = info
# else
  call dposv(uplo,n_,nrhs_,a,lda_,b,ldb_,info_)
# endif
end subroutine dposv_

function dlamch_(cmach)
  use Definitions, only: wp, BLASR8
  implicit none
  real(kind=wp) :: dlamch_
  character :: cmach
  real(kind=BLASR8), external :: dlamch
  dlamch_ = dlamch(cmach)
end function

function ilaenv_(ispec_,nm,opts,n1_,n2_,n3_,n4_)
  use Definitions, only: iwp, BLASInt
  implicit none
  integer(kind=iwp) :: ilaenv_
  integer(kind=iwp) :: ispec_, n1_, n2_, n3_, n4_
  character(len=*) :: nm, opts
  integer(kind=BLASInt), external :: ilaenv
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: ispec, n1, n2, n3, n4
  ispec = int(ispec_,kind=BLASInt)
  n1 = int(n1_,kind=BLASInt)
  n2 = int(n2_,kind=BLASInt)
  n3 = int(n3_,kind=BLASInt)
  n4 = int(n4_,kind=BLASInt)
  ilaenv_ = int(ilaenv(ispec,nm,opts,n1,n2,n3,n4),kind=BLASInt)
# else
  ilaenv_ = ilaenv(ispec_,nm,opts,n1_,n2_,n3_,n4_)
# endif
end function ilaenv_

subroutine dsterf_(n_,d,e,info_)
  use Definitions, only: wp, iwp
  _BLAS_INT_use_
  implicit none
  integer(kind=iwp) :: n_, info_
  real(kind=wp) :: d(*), e(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: info, n
  n = int(n_,kind=BLASInt)
  call dsterf(n,d,e,info)
  info_ = info
# else
  call dsterf(n_,d,e,info_)
# endif
end subroutine dsterf_

subroutine dsteqr_(compz,n_,d,e,z,ldz_,work,info_)
  use Definitions, only: wp, iwp
  _BLAS_INT_use_
  implicit none
  character :: compz
  integer(kind=iwp) :: n_, ldz_, info_
  real(kind=wp) :: d(*), e(*), z(ldz_,*), work(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: info, ldz, n
  n = int(n_,kind=BLASInt)
  ldz = int(ldz_,kind=BLASInt)
  call dsteqr(compz,n,d,e,z,ldz,work,info)
  info_ = info
# else
  call dsteqr(compz,n_,d,e,z,ldz_,work,info_)
# endif
end subroutine dsteqr_

subroutine dlascl_(tp,kl_,ku_,cfrom,cto,m_,n_,a,lda_,info_)
  use Definitions, only: wp, iwp
  _BLAS_INT_use_
  implicit none
  character :: tp
  integer(kind=iwp) :: kl_, ku_, m_, n_, lda_, info_
  real(kind=wp) :: cfrom, cto, a(lda_,*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: info, kl, ku, lda, m, n
  kl = int(kl_,kind=BLASInt)
  ku = int(ku_,kind=BLASInt)
  m = int(m_,kind=BLASInt)
  n = int(n_,kind=BLASInt)
  lda = int(lda_,kind=BLASInt)
  call dlascl(tp,kl,ku,cfrom,cto,m,n,a,lda,info)
  info_ = info
# else
  call dlascl(tp,kl_,ku_,cfrom,cto,m_,n_,a,lda_,info_)
# endif
end subroutine dlascl_

subroutine dorgtr_(uplo,n_,a,lda_,tau,work,lwork_,info_)
  use Definitions, only: wp, iwp
  _BLAS_INT_use_
  implicit none
  character :: uplo
  integer(kind=iwp) :: n_, lda_, lwork_, info_
  real(kind=wp) :: a(lda_,*), tau(*), work(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: info, lda, lwork, n
  n = int(n_,kind=BLASInt)
  lda = int(lda_,kind=BLASInt)
  lwork = int(lwork_,kind=BLASInt)
  call dorgtr(uplo,n,a,lda,tau,work,lwork,info)
  info_ = info
# else
  call dorgtr(uplo,n_,a,lda_,tau,work,lwork_,info_)
# endif
end subroutine dorgtr_

subroutine zgesvd_(jobu,jobvt,m_,n_,a,lda_,s,u,ldu_,vt,ldvt_,work,lwork_,rwork,info_)
  use Definitions, only: wp, iwp
  _BLAS_INT_use_
  implicit none
  character :: jobu, jobvt
  integer(kind=iwp) :: m_, n_, lda_, ldu_, ldvt_, lwork_, info_
  complex(kind=wp) :: a(lda_,*), u(ldu_,*), vt(ldvt_,*), work(*)
  real(kind=wp) :: s(*), rwork(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: info, lda, ldu, ldvt, lwork, m, n
  m = int(m_,kind=BLASInt)
  n = int(n_,kind=BLASInt)
  lda = int(lda_,kind=BLASInt)
  ldu = int(ldu_,kind=BLASInt)
  ldvt = int(ldvt_,kind=BLASInt)
  lwork = int(lwork_,kind=BLASInt)
  call zgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,rwork,info)
  info_ = info
# else
  call zgesvd(jobu,jobvt,m_,n_,a,lda_,s,u,ldu_,vt,ldvt_,work,lwork_,rwork,info_)
# endif
end subroutine zgesvd_

