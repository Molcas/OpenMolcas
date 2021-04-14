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

subroutine DSYEV_VV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO)
!vv identical to DSYEV. workaround for casvb code, which can't stand MKL
!
!  -- LAPACK driver routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999

!     .. Scalar Arguments ..
character JOBZ, UPLO
integer INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
real*8 A(LDA,*), W(*), WORK(*)
!     ..
!
!  Purpose
!  =======
!
!  DSYEV computes all eigenvalues and, optionally, eigenvectors of a
!  real symmetric matrix A.
!
!  Arguments
!  =========
!
!  JOBZ    (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only;
!          = 'V':  Compute eigenvalues and eigenvectors.
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) REAL*8 array, dimension (LDA, N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the
!          leading N-by-N upper triangular part of A contains the
!          upper triangular part of the matrix A.  If UPLO = 'L',
!          the leading N-by-N lower triangular part of A contains
!          the lower triangular part of the matrix A.
!          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
!          orthonormal eigenvectors of the matrix A.
!          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
!          or the upper triangle (if UPLO='U') of A, including the
!          diagonal, is destroyed.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  W       (output) REAL*8 array, dimension (N)
!          If INFO = 0, the eigenvalues in ascending order.
!
!  WORK    (workspace/output) REAL*8 array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The length of the array WORK.  LWORK >= max(1,3*N-1).
!          For optimal efficiency, LWORK >= (NB+2)*N,
!          where NB is the blocksize for DSYTRD returned by ILAENV.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the algorithm failed to converge; i
!                off-diagonal elements of an intermediate tridiagonal
!                form did not converge to zero.
!
!  =====================================================================
!
!     .. Parameters ..
real*8 ZERO, ONE
parameter(ZERO=0.0d0,ONE=1.0d0)
!     ..
!     .. Local Scalars ..
logical LOWER, LQUERY, WANTZ
integer IINFO, IMAX, INDE, INDTAU, INDWRK, ISCALE, LLWORK, LWKOPT, NB
real*8 ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM
!     ..
!     .. External Functions ..
logical LSAME
integer ILAENV_
real*8 DLAMCH_, DLANSY_
external LSAME, ILAENV_, DLAMCH_, DLANSY_
!     ..
!     .. External Subroutines ..
external DLASCL_, DORGTR_, DSCAL_, DSTEQR_, DSTERF_, DSYTRD_, XERBLA
!     ..
!     .. Intrinsic Functions ..
intrinsic MAX, SQRT
!     ..
!     .. Executable Statements ..

! Test the input parameters.

!vv make compiler happy
LWKOPT = -9999999
!vv
WANTZ = LSAME(JOBZ,'V')
LOWER = LSAME(UPLO,'L')
LQUERY = (LWORK == -1)

INFO = 0
if (.not.(WANTZ .or. LSAME(JOBZ,'N'))) then
  INFO = -1
else if (.not.(LOWER .or. LSAME(UPLO,'U'))) then
  INFO = -2
else if (N < 0) then
  INFO = -3
else if (LDA < max(1,N)) then
  INFO = -5
else if (LWORK < max(1,3*N-1) .and. .not. LQUERY) then
  INFO = -8
end if

if (INFO == 0) then
  NB = ILAENV_(1,'DSYTRD',UPLO,N,-1,-1,-1)
  LWKOPT = max(1,(NB+2)*N)
  WORK(1) = LWKOPT
end if

if (INFO /= 0) then
  call XERBLA('DSYEV ',-INFO)
  return
else if (LQUERY) then
  return
end if

! Quick return if possible

if (N == 0) then
  WORK(1) = 1
  return
end if

if (N == 1) then
  W(1) = A(1,1)
  WORK(1) = 3
  if (WANTZ) A(1,1) = ONE
  return
end if

! Get machine constants.

SAFMIN = DLAMCH_('Safe minimum')
EPS = DLAMCH_('Precision')
SMLNUM = SAFMIN/EPS
BIGNUM = ONE/SMLNUM
RMIN = sqrt(SMLNUM)
RMAX = sqrt(BIGNUM)

! Scale matrix to allowable range, if necessary.

ANRM = DLANSY_('M',UPLO,N,A,LDA,WORK)
ISCALE = 0
if (ANRM > ZERO .and. ANRM < RMIN) then
  ISCALE = 1
  SIGMA = RMIN/ANRM
else if (ANRM > RMAX) then
  ISCALE = 1
  SIGMA = RMAX/ANRM
end if
if (ISCALE == 1) call dlascl_(UPLO,0,0,ONE,SIGMA,N,N,A,LDA,INFO)

! Call DSYTRD to reduce symmetric matrix to tridiagonal form.

INDE = 1
INDTAU = INDE+N
INDWRK = INDTAU+N
LLWORK = LWORK-INDWRK+1
call dsytrd_(UPLO,N,A,LDA,W,WORK(INDE),WORK(INDTAU),WORK(INDWRK),LLWORK,IINFO)

! For eigenvalues only, call DSTERF.  For eigenvectors, first call
! DORGTR to generate the orthogonal matrix, then call DSTEQR.

if (.not. WANTZ) then
  call dsterf_(N,W,WORK(INDE),INFO)
else
  call dorgtr_(UPLO,N,A,LDA,WORK(INDTAU),WORK(INDWRK),LLWORK,IINFO)
  call dsteqr_(JOBZ,N,W,WORK(INDE),A,LDA,WORK(INDTAU),INFO)
end if

! If matrix was scaled, then rescale eigenvalues appropriately.

if (ISCALE == 1) then
  if (INFO == 0) then
    IMAX = N
  else
    IMAX = INFO-1
  end if
  call DSCAL_(IMAX,ONE/SIGMA,W,1)
end if

! Set WORK(1) to optimal workspace size.

WORK(1) = LWKOPT

return

! End of DSYEV_VV

end subroutine DSYEV_VV
