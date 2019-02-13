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
! dlansy_
! dsterf_
! dsteqr_
! dlascl_
! dorgtr_

! For procedures known to raise floating point exceptions in the test suite,
! disable exception trapping locally: three pieces of code are needed
#ifdef _FPE_TRAP_
#  define _FPE_TRAP_use_ \
  use, intrinsic :: IEEE_Exceptions
#  define _FPE_TRAP_init_ \
  type(IEEE_Status_Type) :: IEEE_Status ;\
  call IEEE_Get_Status(IEEE_Status) ;\
  call IEEE_Set_Halting_Mode(IEEE_Usual,.false._4)
#  define _FPE_TRAP_end_ \
  call IEEE_Set_Status(IEEE_Status)
#else
#  define _FPE_TRAP_use_ !
#  define _FPE_TRAP_init_ !
#  define _FPE_TRAP_end_ !
#endif

! Set the appropriate integer size of the library interface and specify
! if integer conversion will be needed.
#ifdef LINALG_I4
#  define LAPACKINT INTEGER*4
#  ifdef _I8_
#    define MOLCAS_TO_LAPACK_INT
#  endif
#else
#  define LAPACKINT INTEGER*8
#endif

subroutine dopmtr_(side,uplo,trans,m_,n_,ap,tau,c,ldc_,work,info_)
  implicit none
  character          side, trans, uplo
  integer            info_, ldc_, m_, n_
  real*8             ap( * ), c( ldc_, * ), tau( * ), work( * )
#ifdef MOLCAS_TO_LAPACK_INT
  LAPACKINT          info, ldc, m, n
  m=m_
  n=n_
  ldc=ldc_
  call dopmtr(side,uplo,trans,m,n,ap,tau,c,ldc,work,info)
  info_=info
#else
  call dopmtr(side,uplo,trans,m_,n_,ap,tau,c,ldc_,work,info_)
#endif
end subroutine

subroutine dspgv_(itype_,jobz,uplo,n_,ap,bp,w,z,ldz_,work,info_)
  implicit none
  character          jobz, uplo
  integer            info_, itype_, ldz_, n_
  real*8             ap( * ), bp( * ), w( * ), work( * ), z( ldz_, * )
#ifdef MOLCAS_TO_LAPACK_INT
  LAPACKINT          info, itype, ldz, n
  itype=itype_
  n=n_
  ldz=ldz_
  call dspgv(itype,jobz,uplo,n,ap,bp,w,z,ldz,work,info)
  info_=info
#else
  call dspgv(itype_,jobz,uplo,n_,ap,bp,w,z,ldz_,work,info_)
#endif
end subroutine

subroutine dsptrd_( uplo, n_, ap, d, e, tau, info_ )
  implicit none
  character          uplo
  integer            info_, n_
  real*8             ap( * ), d( * ), e( * ), tau( * )
#ifdef MOLCAS_TO_LAPACK_INT
  LAPACKINT          info, n
  n=n_
  call dsptrd( uplo, n, ap, d, e, tau, info )
  info_=info
#else
  call dsptrd( uplo, n_, ap, d, e, tau, info_ )
#endif
end subroutine

subroutine dstevr_(jobz,range,n_,d,e,vl,vu,il_,iu_,abstol, &
      &            m_,w,z,ldz_,isuppz_,work,lwork_,iwork_, &
      &            liwork_, info_ )
  _FPE_TRAP_use_
  implicit none
  character          jobz, range
  integer            il_, info_, iu_, ldz_, liwork_, lwork_, m_, n_
  real*8             abstol, vl, vu
  integer            isuppz_( * ), iwork_( * )
  real*8             d( * ), e( * ), w( * ), work( * ), z( ldz_, * )
#ifdef MOLCAS_TO_LAPACK_INT
  LAPACKINT          il, info, iu, ldz, liwork, lwork, m, n
  LAPACKINT, allocatable :: isuppz(:), iwork(:)
  integer :: i
  _FPE_TRAP_init_
  n=n_
  il=il_
  iu=iu_
  m=m_
  ldz=ldz_
  liwork=liwork_
  allocate(isuppz(2*max(1,n)))
  allocate(iwork(liwork))
  call dstevr( jobz, range, n, d, e, vl, vu, il, iu, abstol, &
      &        m, w, z, ldz, isuppz, work, lwork, iwork, &
      &        liwork, info )
  deallocate(iwork)
  do i=1,2*max(1,m)
    isuppz_(i)=isuppz(i)
  end do
  deallocate(isuppz)
  info_=info
#else
  _FPE_TRAP_init_
  call dstevr( jobz, range, n_, d, e, vl, vu, il_, iu_, abstol, &
      &        m_, w, z, ldz_, isuppz_, work, lwork_, iwork_, &
      &        liwork_, info_ )
#endif
  _FPE_TRAP_end_
end subroutine

subroutine dgetrs_(trans,n_,nrhs_,a,lda_,ipiv_,b,ldb_,info_)
  implicit none
  character          trans
  integer            info_, lda_, ldb_, n_, nrhs_
  integer            ipiv_( * )
  real*8             a( lda_, * ), b( ldb_, * )
#ifdef MOLCAS_TO_LAPACK_INT
  LAPACKINT          info, lda, ldb, n, nrhs
  LAPACKINT, allocatable :: ipiv(:)
  integer :: i
  n=n_
  nrhs=nrhs_
  lda=lda_
  ldb=ldb_
  allocate(ipiv(n))
  do i=1,n
    ipiv(i)=ipiv_(i)
  end do
  call dgetrs(trans,n,nrhs,a,lda,ipiv,b,ldb,info)
  deallocate(ipiv)
  info_=info
#else
  call dgetrs(trans,n_,nrhs_,a,lda_,ipiv_,b,ldb_,info_)
#endif
end subroutine

subroutine dspev_(jobz,uplo,n_,ap,w,z,ldz_,work,info_)
  implicit none
  character          jobz, uplo
  integer            info_, ldz_, n_
  real*8             ap( * ), w( * ), work( * ), z( ldz_, * )
#ifdef MOLCAS_TO_LAPACK_INT
  LAPACKINT          info, ldz, n
  n=n_
  ldz=ldz_
  call dspev(jobz,uplo,n,ap,w,z,ldz,work,info)
  info_=info
#else
  call dspev(jobz,uplo,n_,ap,w,z,ldz_,work,info_)
#endif
end subroutine

subroutine dgeev_( jobvl, jobvr, n_, a, lda_, wr, wi, vl, ldvl_, &
  &                   vr, ldvr_, work, lwork_, info_ )
  implicit none
  character          jobvl, jobvr
  integer            info_, lda_, ldvl_, ldvr_, lwork_, n_
  real*8             a(lda_,*),wr(*),wi(*),vl(ldvl_,*),vr(ldvr_,*),work(*)
#ifdef MOLCAS_TO_LAPACK_INT
  LAPACKINT          info, lda, ldvl, ldvr, lwork, n
  lda=lda_
  ldvl=ldvl_
  ldvr=ldvr_
  lwork=lwork_
  n=n_
  call dgeev( jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, &
      &                   vr, ldvr, work, lwork, info )
  info_=info
#else
  call dgeev( jobvl, jobvr, n_, a, lda_, wr, wi, vl, ldvl_, &
      &                   vr, ldvr_, work, lwork_, info_ )
#endif
end subroutine

subroutine dgesvd_( jobu, jobvt, m_, n_, a, lda_, s, u, ldu_, &
  &                   vt, ldvt_, work, lwork_, info_ )
  implicit none
  character          jobu, jobvt
  integer            info_, lda_, ldu_, ldvt_, lwork_, m_, n_
  real*8             a(lda_,*),s(*),u(ldu_,*),vt(ldvt_,*),work(*)
#ifdef MOLCAS_TO_LAPACK_INT
  LAPACKINT          info, lda, ldu, ldvt, lwork, m, n
  lda=lda_
  ldu=ldu_
  ldvt=ldvt_
  lwork=lwork_
  m=m_
  n=n_
  call dgesvd( jobu, jobvt, m, n, a, lda, s, u, ldu, &
      &                   vt, ldvt, work, lwork, info )
  info_=info
#else
  call dgesvd( jobu, jobvt, m_, n_, a, lda_, s, u, ldu_, &
      &                   vt, ldvt_, work, lwork_, info_ )
#endif
end subroutine

subroutine dgetrf_( m_,n_,a,lda_,ipiv_,info_ )
  implicit none
  integer            info_,lda_,m_,n_
  integer            ipiv_( * )
  real*8             a( lda_, * )
#ifdef MOLCAS_TO_LAPACK_INT
  LAPACKINT          info, lda, m, n
  LAPACKINT, allocatable :: ipiv(:)
  integer :: i
  allocate(ipiv(n))
  do i=1,n
    ipiv(i)=ipiv_(i)
  end do
  call dgetrf( m,n,a,lda,ipiv,info )
  deallocate(ipiv)
  info_=info
#else
  call dgetrf( m_,n_,a,lda_,ipiv_,info_ )
#endif
end subroutine

subroutine dgesv_( n_, nrhs_, a, lda_, ipiv_, b, ldb_, info_ )
  implicit none
  integer            info_, lda_, ldb_, n_, nrhs_
  integer            ipiv_( * )
  real*8             a( lda_, * ), b( ldb_, * )
#ifdef MOLCAS_TO_LAPACK_INT
  LAPACKINT          info, lda, ldb, n, nrhs
  LAPACKINT, allocatable :: ipiv(:)
  integer :: i
  n=n_
  nrhs=nrhs_
  lda=lda_
  ldb=ldb_
  allocate(ipiv(n))
  do i=1,n
    ipiv(i)=ipiv_(i)
  end do
  call dgesv( n, nrhs, a, lda, ipiv, b, ldb, info )
  deallocate(ipiv)
#else
  call dgesv( n_, nrhs_, a, lda_, ipiv_, b, ldb_, info_ )
#endif
end subroutine

subroutine dsyev_( jobz, uplo, n_, a, lda_, w, work, lwork_, info_ )
  implicit none
  character          jobz, uplo
  integer            info_, lda_, lwork_, n_
  real*8             a( lda_, * ), w( * ), work( * )
#ifdef MOLCAS_TO_LAPACK_INT
  LAPACKINT          info, lda, lwork, n
  n=n_
  lda=lda_
  lwork=lwork_
  call dsyev( jobz, uplo, n, a, lda, w, work, lwork, info )
  info_=info
#else
  call dsyev( jobz, uplo, n_, a, lda_, w, work, lwork_, info_ )
#endif
end subroutine

subroutine dsyevr_( jobz, range, uplo, n_, a, lda_, vl, vu, il_, iu_, &
      &            abstol, m_, w, z, ldz_, isuppz_, work, lwork_, &
      &            iwork_, liwork_, info_ )
  implicit none
  character          jobz, range, uplo
  integer            il_, info_, iu_, lda_, ldz_, liwork_, lwork_, m_, n_
  real*8             abstol, vl, vu
  integer            isuppz_( * ), iwork_( * )
  real*8             a( lda_, * ), w( * ), work( * ), z( ldz_, * )
#ifdef MOLCAS_TO_LAPACK_INT
  LAPACKINT          il, info, iu, lda, ldz, liwork, lwork, m, n
  LAPACKINT, allocatable :: isuppz(:), iwork(:)
  integer :: i
  n=n_
  lda=lda_
  il=il_
  iu=iu_
  m=m_
  ldz=ldz_
  lwork=lwork_
  liwork=liwork_
  allocate(isuppz(2*max(1,n)))
  allocate(iwork(max(1,liwork)))
  call dsyevr( jobz, range, uplo, n, a, lda, vl, vu, il, iu, &
      &        abstol, m, w, z, ldz, isuppz, work, lwork, &
      &        iwork, liwork, info )
  do i=1,2*max(1,m)
    isuppz_(i)=isuppz(i)
  end do
  deallocate(isuppz)
  iwork_(1)=iwork(1)
  deallocate(iwork)
  info_=info
#else
  call dsyevr( jobz, range, uplo, n_, a, lda_, vl, vu, il_, iu_, &
      &        abstol, m_, w, z, ldz_, isuppz_, work, lwork_, &
      &        iwork_, liwork_, info_ )
#endif
end subroutine

subroutine dgetri_( n_, a, lda_, ipiv_, work, lwork_, info_ )
  implicit none
  integer            info_, lda_, lwork_, n_
  integer            ipiv_( * )
  real*8             a( lda_, * ), work( * )
#ifdef MOLCAS_TO_LAPACK_INT
  LAPACKINT            info, lda, lwork, n
  LAPACKINT, allocatable :: ipiv(:)
  integer :: i
  n=n_
  lda=lda_
  lwork=lwork_
  allocate(ipiv(n))
  do i=1,n
    ipiv(i)=ipiv_(i)
  end do
  call dgetri( n, a, lda, ipiv, work, lwork, info )
  deallocate(ipiv)
  info_=info
#else
  call dgetri( n_, a, lda_, ipiv_, work, lwork_, info_ )
#endif
end subroutine

subroutine zhpev_( jobz, uplo, n_, ap, w, z, ldz_, work, rwork, info_ )
  implicit none
  character          jobz, uplo
  integer            info_, ldz_, n_
  real*8             rwork( * ), w( * )
  complex*16         ap( * ), work( * ), z( ldz_, * )
#ifdef MOLCAS_TO_LAPACK_INT
  LAPACKINT          info, ldz, n
  n=n_
  ldz=ldz_
  call zhpev( jobz, uplo, n, ap, w, z, ldz, work, rwork, info )
  info_=info
#else
  call zhpev( jobz, uplo, n_, ap, w, z, ldz_, work, rwork, info_ )
#endif
end subroutine

subroutine dsygv_( itype_, jobz, uplo, n_, a, lda_, b, ldb_, w, &
      &            work, lwork_, info_ )
  implicit none
  character          jobz, uplo
  integer            info_, itype_, lda_, ldb_, lwork_, n_
  real*8             a( lda_, * ), b( ldb_, * ), w( * ), work( * )
#ifdef MOLCAS_TO_LAPACK_INT
  LAPACKINT          info, itype, lda, ldb, lwork, n
  itype=itype_
  n=n_
  lda=lda_
  ldb=ldb_
  lwork=lwork_
  call dsygv( itype, jobz, uplo, n, a, lda, b, ldb, w, work, &
      &       lwork, info )
  info_=info
#else
  call dsygv( itype_, jobz, uplo, n_, a, lda_, b, ldb_, w, work, &
      &       lwork_, info_ )
#endif
end subroutine

subroutine dpotrf_( uplo, n_, a, lda_, info_ )
  implicit none
  character          uplo
  integer            info_, lda_, n_
  real*8             a( lda_, * )
#ifdef MOLCAS_TO_LAPACK_INT
  LAPACKINT          info, lda, n
  n=n_
  lda=lda_
  call dpotrf( uplo, n, a, lda, info )
  info_=info
#else
  call dpotrf( uplo, n_, a, lda_, info_ )
#endif
end subroutine

subroutine dgels_( trans, m_, n_, nrhs_, a, lda_, b, ldb_, work, lwork_, info_ )
  implicit none
  character          trans
  integer            info_, lda_, ldb_, lwork_, m_, n_, nrhs_
  real*8             a( lda_, * ), b( ldb_, * ), work( * )
#ifdef MOLCAS_TO_LAPACK_INT
  LAPACKINT          info, lda, ldb, lwork, m, n, nrhs
  m=m_
  n=n_
  nrhs=nrhs_
  lda=lda_
  ldb=ldb_
  lwork=lwork_
  call dgels( trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info )
  info_=info
#else
  call dgels( trans, m_, n_, nrhs_, a, lda_, b, ldb_, work, lwork_, info_ )
#endif
end subroutine

subroutine dsytrd_( uplo, n_, a, lda_, d, e, tau, work, lwork_, info_ )
  implicit none
  character          uplo
  integer            info_, lda_, lwork_, n_
  real*8             a( lda_, * ), d( * ), e( * ), tau( * ), work (*)
#ifdef MOLCAS_TO_LAPACK_INT
  LAPACKINT          info, lda, lwork, n
  n=n_
  lda=lda_
  lwork=lwork_
  call dsytrd( uplo, n, a, lda, d, e, tau, work, lwork, info )
  info_=info
#else
  call dsytrd( uplo, n_, a, lda_, d, e, tau, work, lwork_, info_ )
#endif
end subroutine

subroutine dposv_( uplo, n_, nrhs_, a, lda_, b, ldb_, info_ )
  implicit none
  character          uplo
  integer            info_, lda_, ldb_, n_, nrhs_
  real*8             a( lda_, * ), b( ldb_, * )
#ifdef MOLCAS_TO_LAPACK_INT
  LAPACKINT          info, lda, ldb, n, nrhs
  n=n_
  nrhs=nrhs_
  lda=lda_
  ldb=ldb_
  call dposv( uplo, n, nrhs, a, lda, b, ldb, info )
  info_=info
#else
  call dposv( uplo, n_, nrhs_, a, lda_, b, ldb_, info_ )
#endif
end subroutine

real*8 function dlamch_( cmach )
  implicit none
  character :: cmach
  real*8, external :: dlamch
  dlamch_ = dlamch(cmach)
end function

integer function ilaenv_( ispec_, name, opts, n1_, n2_, n3_, n4_ )
  implicit none
  character(*)    name, opts
  integer            ispec_, n1_, n2_, n3_, n4_
#ifdef MOLCAS_TO_LAPACK_INT
  LAPACKINT          ispec, n1, n2, n3, n4
  LAPACKINT, external :: ilaenv
  ispec=ispec_
  n1=n1_
  n2=n2_
  n3=n3_
  n4=n4_
  ilaenv_ = ilaenv( ispec, name, opts, n1, n2, n3, n4 )
#else
  integer, external :: ilaenv
  ilaenv_ = ilaenv( ispec_, name, opts, n1_, n2_, n3_, n4_ )
#endif
end function

real*8 function dlansy_( norm, uplo, n_, a, lda_, work )
  implicit none
  character          norm, uplo
  integer            lda_, n_
  real*8             a( lda_, * ), work( * )
  real*8, external :: dlansy
#ifdef MOLCAS_TO_LAPACK_INT
  LAPACKINT          lda, n
  n=n_
  lda=lda_
  dlansy_ = dlansy( norm, uplo, n, a, lda, work )
#else
  dlansy_ = dlansy( norm, uplo, n_, a, lda_, work )
#endif
end function

subroutine dsterf_( n_, d, e, info_ )
  implicit none
  integer            info_, n_
  real*8             d( * ), e( * )
#ifdef MOLCAS_TO_LAPACK_INT
  LAPACKINT          info, n
  n=n_
  call dsterf( n, d, e, info )
  info_=info
#else
  call dsterf( n_, d, e, info_ )
#endif
end subroutine

subroutine dsteqr_( compz, n_, d, e, z, ldz_, work, info_ )
  implicit none
  character          compz
  integer            info_, ldz_, n_
  real*8             d( * ), e( * ), work( * ), z( ldz_, * )
#ifdef MOLCAS_TO_LAPACK_INT
  LAPACKINT          info, ldz, n
  n=n_
  ldz=ldz_
  call dsteqr( compz, n, d, e, z, ldz, work, info )
  info_=info
#else
  call dsteqr( compz, n_, d, e, z, ldz_, work, info_ )
#endif
end subroutine

subroutine dlascl_( type, kl_, ku_, cfrom, cto, m_, n_, a, lda_, info_ )
  implicit none
  character          type
  integer            info_, kl_, ku_, lda_, m_, n_
  real*8             cfrom, cto
  real*8             a( lda_, * )
#ifdef MOLCAS_TO_LAPACK_INT
  LAPACKINT          info, kl, ku, lda, m, n
  kl=kl_
  ku=ku_
  m=m_
  n=n_
  lda=lda_
  call dlascl( type, kl, ku, cfrom, cto, m, n, a, lda, info )
  info_=info
#else
  call dlascl( type, kl_, ku_, cfrom, cto, m_, n_, a, lda_, info_ )
#endif
end subroutine

subroutine dorgtr_( uplo, n_, a, lda_, tau, work, lwork_, info_ )
  implicit none
  character          uplo
  integer            info_, lda_, lwork_, n_
  real*8             a( lda_, * ), tau( * ), work( * )
#ifdef MOLCAS_TO_LAPACK_INT
  LAPACKINT          info, lda, lwork, n
  n=n_
  lda=lda_
  lwork=lwork_
  call dorgtr( uplo, n, a, lda, tau, work, lwork, info )
  info_=info
#else
  call dorgtr( uplo, n_, a, lda_, tau, work, lwork_, info_ )
#endif
end subroutine
