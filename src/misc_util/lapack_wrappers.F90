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
! dgecon_
! dgees_
! dgeev_
! dgels_
! dgesv_
! dgesvd_
! dgetrf_
! dgetri_
! dgetrs_
! dlange_
! dlascl_
! dopmtr_
! dorgtr_
! dposv_
! dpotrf_
! dspev_
! dspgv_
! dsptrd_
! dsteqr_
! dsterf_
! dstevr_
! dsyev_
! dsyevr_
! dsygv_
! dsytrd_
! ilaenv_
! zgesvd_
! zhpev_

! Specify if integer (and logical) conversion will be needed.
! (real conversion is not implemented yet)
#if defined(LINALG_I4) && defined(_I8_)
# define MOLCAS_TO_BLAS_INT
# define _BLAS_INT_use_ use Definitions, only: BLASInt
# define _BLAS_INT_stdalloc_ \
  use stdalloc, only: mma_allocate, mma_deallocate
#else
# define _BLAS_INT_use_
# define _BLAS_INT_stdalloc_ !
#endif

#include "intent.fh"

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

subroutine dgecon_(norm,n_,a,lda_,anorm,rcond,work,iwork_,info_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  _BLAS_INT_stdalloc_
  implicit none
  character, intent(in) :: norm
  integer(kind=iwp), intent(in) :: n_, lda_
  real(kind=BLASR8), intent(in) :: a(lda_,*), anorm
  real(kind=BLASR8), intent(out) :: rcond
  real(kind=BLASR8), intent(_OUT_) :: work(*)
  integer(kind=iwp), intent(_OUT_) :: iwork_(*)
  integer(kind=iwp), intent(out) :: info_
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: info, lda, n
  integer(kind=BLASInt), allocatable :: iwork(:)
  n = int(n_,kind=BLASInt)
  lda = int(lda_,kind=BLASInt)
  call mma_allocate(iwork,n_,label='iwork')
  call dgecon(norm,n,a,lda,anorm,rcond,work,iwork,info)
  call mma_deallocate(iwork)
  info_ = info
# else
  call dgecon(norm,n_,a,lda_,anorm,rcond,work,iwork_,info_)
# endif
end subroutine dgecon_

subroutine dgees_(jobvs,sort,slct_,n_,a,lda_,sdim_,wr,wi,vs,ldvs_,work,lwork_,bwork_,info_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  _BLAS_INT_stdalloc_
  implicit none
  character, intent(in) :: jobvs, sort
  integer(kind=iwp), intent(in) :: n_, lda_, ldvs_, lwork_
  real(kind=BLASR8), intent(inout) :: a(lda_,*)
  integer(kind=iwp), intent(out) :: sdim_, info_
  real(kind=BLASR8), intent(_OUT_) :: wr(*), wi(*), vs(ldvs_,*), work(*)
  logical(kind=iwp), external :: slct_
  logical(kind=iwp), intent(_OUT_) :: bwork_(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: n, lda, sdim, ldvs, lwork, info
  logical(kind=BLASInt), allocatable :: bwork(:)
  n = int(n_,kind=BLASInt)
  lda = int(lda_,kind=BLASInt)
  ldvs = int(ldvs_,kind=BLASInt)
  lwork = int(lwork_,kind=BLASInt)
  call mma_allocate(bwork,n_,label='bwork')
  call dgees(jobvs,sort,slct,n,a,lda,sdim,wr,wi,vs,ldvs,work,lwork,bwork,info)
  bwork_(1) = bwork(1)
  call mma_deallocate(bwork)
  sdim_ = sdim
  info_ = info
  contains
  function slct(a,b)
  logical(kind=BLASInt) :: slct
  real(kind=BLASR8), intent(in) :: a, b
  slct = logical(slct_(a,b),kind=BLASInt)
  end function slct
# else
  call dgees(jobvs,sort,slct_,n_,a,lda_,sdim_,wr,wi,vs,ldvs_,work,lwork_,bwork_,info_)
# endif
end subroutine dgees_

subroutine dgeev_(jobvl,jobvr,n_,a,lda_,wr,wi,vl,ldvl_,vr,ldvr_,work,lwork_,info_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  character, intent(in) :: jobvl, jobvr
  integer(kind=iwp), intent(in) :: n_, lda_, ldvl_, ldvr_, lwork_
  real(kind=BLASR8), intent(inout) :: a(lda_,*)
  real(kind=BLASR8), intent(_OUT_) :: wr(*), wi(*), vl(ldvl_,*), vr(ldvr_,*), work(*)
  integer(kind=iwp), intent(out) :: info_
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

subroutine dgels_(trans,m_,n_,nrhs_,a,lda_,b,ldb_,work,lwork_,info_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  character, intent(in) :: trans
  integer(kind=iwp), intent(in) :: m_, n_, nrhs_, lda_, ldb_, lwork_
  real(kind=BLASR8), intent(inout) :: a(lda_,*), b(ldb_,*)
  real(kind=BLASR8), intent(_OUT_) :: work(*)
  integer(kind=iwp), intent(out) :: info_
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

subroutine dgesv_(n_,nrhs_,a,lda_,ipiv_,b,ldb_,info_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  _BLAS_INT_stdalloc_
  implicit none
  integer(kind=iwp), intent(in) :: n_, nrhs_, lda_, ldb_
  real(kind=BLASR8), intent(inout) :: a(lda_,*), b(ldb_,*)
  integer(kind=iwp), intent(_OUT_) :: ipiv_(*)
  integer(kind=iwp), intent(out) :: info_
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: info, lda, ldb, n, nrhs
  integer(kind=BLASInt), allocatable :: ipiv(:)
  n = int(n_,kind=BLASInt)
  nrhs = int(nrhs_,kind=BLASInt)
  lda = int(lda_,kind=BLASInt)
  ldb = int(ldb_,kind=BLASInt)
  call mma_allocate(ipiv,n_,label='ipiv')
  call dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
  ipiv_(1:n_) = ipiv
  call mma_deallocate(ipiv)
  info_ = info
# else
  call dgesv(n_,nrhs_,a,lda_,ipiv_,b,ldb_,info_)
# endif
end subroutine dgesv_

subroutine dgesvd_(jobu,jobvt,m_,n_,a,lda_,s,u,ldu_,vt,ldvt_,work,lwork_,info_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  character, intent(in) :: jobu, jobvt
  integer(kind=iwp), intent(in) :: m_, n_, lda_, ldu_, ldvt_, lwork_
  real(kind=BLASR8), intent(inout) :: a(lda_,*)
  real(kind=BLASR8), intent(_OUT_) :: s(*), u(ldu_,*), vt(ldvt_,*), work(*)
  integer(kind=iwp), intent(out) :: info_
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
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  _BLAS_INT_stdalloc_
  implicit none
  integer(kind=iwp), intent(in) :: m_, n_, lda_
  real(kind=BLASR8), intent(inout) :: a(lda_,*)
  integer(kind=iwp), intent(_OUT_) :: ipiv_(*)
  integer(kind=iwp), intent(out) :: info_
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=iwp) :: npiv
  integer(kind=BLASInt) :: info, lda, m, n
  integer(kind=BLASInt), allocatable :: ipiv(:)
  m = int(m_,kind=BLASInt)
  n = int(n_,kind=BLASInt)
  lda = int(lda_,kind=BLASInt)
  npiv = min(m_,n_)
  call mma_allocate(ipiv,npiv,label='ipiv')
  call dgetrf(m,n,a,lda,ipiv,info)
  ipiv_(1:npiv) = ipiv
  call mma_deallocate(ipiv)
  info_ = info
# else
  call dgetrf(m_,n_,a,lda_,ipiv_,info_)
# endif
end subroutine dgetrf_

subroutine dgetri_(n_,a,lda_,ipiv_,work,lwork_,info_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  _BLAS_INT_stdalloc_
  implicit none
  integer(kind=iwp), intent(in) :: n_, lda_, ipiv_(*), lwork_
  real(kind=BLASR8), intent(inout) :: a(lda_,*)
  real(kind=BLASR8), intent(_OUT_) :: work(*)
  integer(kind=iwp), intent(out) :: info_
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: info, lda, lwork, n
  integer(kind=BLASInt), allocatable :: ipiv(:)
  n = int(n_,kind=BLASInt)
  lda = int(lda_,kind=BLASInt)
  call mma_allocate(ipiv,n_,label='ipiv')
  ipiv(:) = int(ipiv_(1:n_),kind=BLASInt)
  lwork = int(lwork_,kind=BLASInt)
  call dgetri(n,a,lda,ipiv,work,lwork,info)
  call mma_deallocate(ipiv)
  info_ = info
# else
  call dgetri(n_,a,lda_,ipiv_,work,lwork_,info_)
# endif
end subroutine dgetri_

subroutine dgetrs_(trans,n_,nrhs_,a,lda_,ipiv_,b,ldb_,info_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  _BLAS_INT_stdalloc_
  implicit none
  character :: trans
  integer(kind=iwp), intent(in) :: n_, nrhs_, lda_, ipiv_(*), ldb_
  real(kind=BLASR8), intent(in) :: a(lda_,*)
  real(kind=BLASR8), intent(inout) :: b(ldb_,*)
  integer(kind=iwp), intent(out) :: info_
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: info, lda, ldb, n, nrhs
  integer(kind=BLASInt), allocatable :: ipiv(:)
  n = int(n_,kind=BLASInt)
  nrhs = int(nrhs_,kind=BLASInt)
  lda = int(lda_,kind=BLASInt)
  call mma_allocate(ipiv,n_,label='ipiv')
  ipiv(:) = int(ipiv_(1:n_),kind=BLASInt)
  ldb = int(ldb_,kind=BLASInt)
  call dgetrs(trans,n,nrhs,a,lda,ipiv,b,ldb,info)
  call mma_deallocate(ipiv)
  info_ = info
# else
  call dgetrs(trans,n_,nrhs_,a,lda_,ipiv_,b,ldb_,info_)
# endif
end subroutine dgetrs_

function dlange_(norm,m_,n_,a,lda_,work)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  real(kind=BLASR8) :: dlange_
  character, intent(in) :: norm
  integer(kind=iwp), intent(in) :: m_, n_, lda_
  real(kind=BLASR8), intent(in) :: a(lda_,*)
  real(kind=BLASR8), intent(_OUT_) :: work(*)
  real(kind=BLASR8), external :: dlange
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: lda, m, n
  m = int(m_,kind=BLASInt)
  n = int(n_,kind=BLASInt)
  lda = int(lda_,kind=BLASInt)
  dlange_ = dlange(norm,m,n,a,lda,work)
# else
  dlange_ = dlange(norm,m_,n_,a,lda_,work)
# endif
end function dlange_

subroutine dlascl_(tp,kl_,ku_,cfrom,cto,m_,n_,a,lda_,info_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  character, intent(in) :: tp
  integer(kind=iwp), intent(in) :: kl_, ku_, m_, n_, lda_
  real(kind=BLASR8), intent(in) :: cfrom, cto
  real(kind=BLASR8), intent(inout) :: a(lda_,*)
  integer(kind=iwp), intent(out) :: info_
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

subroutine dopmtr_(side,uplo,trans,m_,n_,ap,tau,c,ldc_,work,info_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  character, intent(in) :: side, uplo, trans
  integer(kind=iwp), intent(in) :: m_, n_, ldc_
  real(kind=BLASR8), intent(in) :: ap(*), tau(*)
  real(kind=BLASR8), intent(inout) :: c(ldc_,*)
  real(kind=BLASR8), intent(_OUT_) :: work(*)
  integer(kind=iwp), intent(out) :: info_
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

subroutine dorgtr_(uplo,n_,a,lda_,tau,work,lwork_,info_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  character :: uplo
  integer(kind=iwp), intent(in) :: n_, lda_, lwork_
  real(kind=BLASR8), intent(inout) :: a(lda_,*), tau(*)
  real(kind=BLASR8), intent(_OUT_) :: work(*)
  integer(kind=iwp), intent(out) :: info_
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

subroutine dposv_(uplo,n_,nrhs_,a,lda_,b,ldb_,info_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  character :: uplo
  integer(kind=iwp), intent(in) :: n_, nrhs_, lda_, ldb_
  real(kind=BLASR8), intent(inout) :: a(lda_,*), b(ldb_,*)
  integer(kind=iwp), intent(out) :: info_
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

subroutine dpotrf_(uplo,n_,a,lda_,info_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  character :: uplo
  integer(kind=iwp), intent(in) :: n_, lda_
  real(kind=BLASR8), intent(inout) :: a(lda_,*)
  integer(kind=iwp), intent(out) :: info_
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

subroutine dspev_(jobz,uplo,n_,ap,w,z,ldz_,work,info_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  character, intent(in) :: jobz, uplo
  integer(kind=iwp), intent(in) :: n_, ldz_
  real(kind=BLASR8), intent(inout) :: ap(*)
  real(kind=BLASR8), intent(_OUT_) :: w(*), z(ldz_,*), work(*)
  integer(kind=iwp), intent(out) :: info_
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

subroutine dspgv_(itype_,jobz,uplo,n_,ap,bp,w,z,ldz_,work,info_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  integer(kind=iwp), intent(in) :: itype_, n_, ldz_
  character, intent(in) :: jobz, uplo
  real(kind=BLASR8), intent(inout) :: ap(*), bp(*)
  real(kind=BLASR8), intent(_OUT_) :: w(*), z(ldz_,*), work(*)
  integer(kind=iwp), intent(out) :: info_
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
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  character, intent(in) :: uplo
  integer(kind=iwp), intent(in) :: n_
  real(kind=BLASR8), intent(inout) :: ap(*)
  real(kind=BLASR8), intent(_OUT_) :: d(*), e(*), tau(*)
  integer(kind=iwp), intent(out) :: info_
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: n, info
  n = int(n_,kind=BLASInt)
  call dsptrd(uplo,n,ap,d,e,tau,info)
  info_ = info
# else
  call dsptrd(uplo,n_,ap,d,e,tau,info_)
# endif
end subroutine dsptrd_

subroutine dsteqr_(compz,n_,d,e,z,ldz_,work,info_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  character, intent(in) :: compz
  integer(kind=iwp), intent(in) :: n_, ldz_
  real(kind=BLASR8), intent(inout) :: d(*), e(*), z(ldz_,*)
  real(kind=BLASR8), intent(_OUT_) :: work(*)
  integer(kind=iwp), intent(out) :: info_
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

subroutine dsterf_(n_,d,e,info_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  integer(kind=iwp), intent(in) :: n_
  real(kind=BLASR8), intent(inout) :: d(*), e(*)
  integer(kind=iwp), intent(out) :: info_
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=BLASInt) :: info, n
  n = int(n_,kind=BLASInt)
  call dsterf(n,d,e,info)
  info_ = info
# else
  call dsterf(n_,d,e,info_)
# endif
end subroutine dsterf_

subroutine dstevr_(jobz,rng,n_,d,e,vl,vu,il_,iu_,abstol,m_,w,z,ldz_,isuppz_,work,lwork_,iwork_,liwork_,info_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  _BLAS_INT_stdalloc_
  _FPE_TRAP_use_
  implicit none
  character, intent(in) :: jobz, rng
  integer(kind=iwp), intent(in) :: n_, il_, iu_, ldz_, lwork_, liwork_
  real(kind=BLASR8), intent(inout) :: d(*), e(*)
  real(kind=BLASR8), intent(in) :: vl, vu, abstol
  integer(kind=iwp), intent(out) :: m_, info_
  real(kind=BLASR8), intent(_OUT_) :: w(*), z(ldz_,*), work(*)
  integer(kind=iwp), intent(_OUT_) :: isuppz_(*), iwork_(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=iwp) :: nsuppz
  integer(kind=BLASInt) :: il, info, iu, ldz, liwork, lwork, m, n
  integer(kind=BLASInt), allocatable :: isuppz(:), iwork(:)
  _FPE_TRAP_init_
  n = int(n_,kind=BLASInt)
  il = int(il_,kind=BLASInt)
  iu = int(iu_,kind=BLASInt)
  ldz = int(ldz_,kind=BLASInt)
  lwork = int(lwork_,kind=BLASInt)
  liwork = int(liwork_,kind=BLASInt)
  nsuppz = 2*max(1,n_)
  call mma_allocate(isuppz,nsuppz,label='isuppz')
  call mma_allocate(iwork,liwork_,label='iwork')
  call dstevr(jobz,rng,n,d,e,vl,vu,il,iu,abstol,m,w,z,ldz,isuppz,work,lwork,iwork,liwork,info)
  m_ = m
  nsuppz = 2*max(1,m_)
  isuppz_(1:nsuppz) = isuppz(1:nsuppz)
  call mma_deallocate(isuppz)
  iwork_(1) = iwork(1)
  call mma_deallocate(iwork)
  info_ = info
# else
  _FPE_TRAP_init_
  call dstevr(jobz,rng,n_,d,e,vl,vu,il_,iu_,abstol,m_,w,z,ldz_,isuppz_,work,lwork_,iwork_,liwork_,info_)
# endif
  _FPE_TRAP_end_
end subroutine dstevr_

subroutine dsyev_(jobz,uplo,n_,a,lda_,w,work,lwork_,info_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  character, intent(in) :: jobz, uplo
  integer(kind=iwp), intent(in) :: n_, lda_, lwork_
  real(kind=BLASR8), intent(inout) :: a(lda_,*)
  real(kind=BLASR8), intent(_OUT_) :: w(*), work(*)
  integer(kind=iwp), intent(out) :: info_
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
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  _BLAS_INT_stdalloc_
  implicit none
  character, intent(in) :: jobz, rng, uplo
  integer(kind=iwp), intent(in) :: n_, lda_, il_, iu_, ldz_, lwork_, liwork_
  real(kind=BLASR8), intent(inout) :: a(lda_,*)
  real(kind=BLASR8), intent(in) :: vl, vu, abstol
  integer(kind=iwp), intent(out) :: m_, info_
  real(kind=BLASR8), intent(_OUT_) :: w(*), z(ldz_,*), work(*)
  integer(kind=iwp), intent(_OUT_) :: isuppz_(*), iwork_(*)
# ifdef MOLCAS_TO_BLAS_INT
  integer(kind=iwp) :: nsuppz
  integer(kind=BLASInt) :: il, info, iu, lda, ldz, liwork, lwork, m, n
  integer(kind=BLASInt), allocatable :: isuppz(:), iwork(:)
  n = int(n_,kind=BLASInt)
  lda = int(lda_,kind=BLASInt)
  il = int(il_,kind=BLASInt)
  iu = int(iu_,kind=BLASInt)
  ldz = int(ldz_,kind=BLASInt)
  lwork = int(lwork_,kind=BLASInt)
  liwork = int(liwork_,kind=BLASInt)
  nsuppz = 2*max(1,n_)
  call mma_allocate(isuppz,nsuppz,label='isuppz')
  call mma_allocate(iwork,max(1,liwork_),label='iwork')
  call dsyevr(jobz,rng,uplo,n,a,lda,vl,vu,il,iu,abstol,m,w,z,ldz,isuppz,work,lwork,iwork,liwork,info)
  m_ = m
  nsuppz = 2*max(1,n_)
  isuppz_(1:nsuppz) = isuppz(1:nsuppz)
  call mma_deallocate(isuppz)
  iwork_(1) = iwork(1)
  call mma_deallocate(iwork)
  info_ = info
# else
  call dsyevr(jobz,rng,uplo,n_,a,lda_,vl,vu,il_,iu_,abstol,m_,w,z,ldz_,isuppz_,work,lwork_,iwork_,liwork_,info_)
# endif
end subroutine dsyevr_

subroutine dsygv_(itype_,jobz,uplo,n_,a,lda_,b,ldb_,w,work,lwork_,info_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  integer(kind=iwp), intent(in) :: itype_, n_, lda_, ldb_, lwork_
  character, intent(in) :: jobz, uplo
  real(kind=BLASR8), intent(inout) :: a(lda_,*), b(ldb_,*)
  real(kind=BLASR8), intent(_OUT_) :: w(*), work(*)
  integer(kind=iwp), intent(out) :: info_
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

subroutine dsytrd_(uplo,n_,a,lda_,d,e,tau,work,lwork_,info_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  character, intent(in) :: uplo
  integer(kind=iwp), intent(in) :: n_, lda_, lwork_
  real(kind=BLASR8), intent(inout) :: a(lda_,*)
  real(kind=BLASR8), intent(_OUT_) :: d(*), e(*), tau(*), work(*)
  integer(kind=iwp), intent(out) :: info_
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

function ilaenv_(ispec_,nm,opts,n1_,n2_,n3_,n4_)
  use Definitions, only: iwp, BLASInt
  implicit none
  integer(kind=iwp) :: ilaenv_
  integer(kind=iwp), intent(in) :: ispec_, n1_, n2_, n3_, n4_
  character(len=*), intent(in) :: nm, opts
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

subroutine zgesvd_(jobu,jobvt,m_,n_,a,lda_,s,u,ldu_,vt,ldvt_,work,lwork_,rwork,info_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  character, intent(in) :: jobu, jobvt
  integer(kind=iwp), intent(in) :: m_, n_, lda_, ldu_, ldvt_, lwork_
  complex(kind=BLASR8), intent(inout) :: a(lda_,*)
  real(kind=BLASR8), intent(_OUT_) :: s(*), rwork(*)
  complex(kind=BLASR8), intent(_OUT_) :: u(ldu_,*), vt(ldvt_,*), work(*)
  integer(kind=iwp), intent(out) :: info_
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

subroutine zhpev_(jobz,uplo,n_,ap,w,z,ldz_,work,rwork,info_)
  use Definitions, only: BLASR8, iwp
  _BLAS_INT_use_
  implicit none
  character, intent(in) :: jobz, uplo
  integer(kind=iwp), intent(in) :: n_, ldz_
  complex(kind=BLASR8), intent(inout) :: ap(*)
  real(kind=BLASR8), intent(_OUT_) :: w(*), rwork(*)
  complex(kind=BLASR8), intent(_OUT_) :: z(ldz_,*), work(*)
  integer(kind=iwp), intent(out) :: info_
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
