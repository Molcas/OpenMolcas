************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE DSYEV_VV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
cvv identical to DSYEV. workaround for casvb code, which can't stand MKL
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      REAL*8   A( LDA, * ), W( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DSYEV computes all eigenvalues and, optionally, eigenvectors of a
*  real symmetric matrix A.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) REAL*8 array, dimension (LDA, N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of A contains the
*          upper triangular part of the matrix A.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of A contains
*          the lower triangular part of the matrix A.
*          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
*          orthonormal eigenvectors of the matrix A.
*          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
*          or the upper triangle (if UPLO='U') of A, including the
*          diagonal, is destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  W       (output) REAL*8 array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*
*  WORK    (workspace/output) REAL*8 array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= max(1,3*N-1).
*          For optimal efficiency, LWORK >= (NB+2)*N,
*          where NB is the blocksize for DSYTRD returned by ILAENV.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the algorithm failed to converge; i
*                off-diagonal elements of an intermediate tridiagonal
*                form did not converge to zero.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL*8   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LOWER, LQUERY, WANTZ
      INTEGER            IINFO, IMAX, INDE, INDTAU, INDWRK, ISCALE,
     $                   LLWORK, LOPT, LWKOPT, NB
      REAL*8   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA,
     $                   SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV_
      REAL*8   DLAMCH_, DLANSY_
      EXTERNAL           LSAME, ILAENV_, DLAMCH_, DLANSY_
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASCL_, DORGTR_, DSCAL_, DSTEQR_, DSTERF_,
     $                   DSYTRD_, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
cvv make compiler happy
      LWKOPT=-9999999
cvv
      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
      LQUERY = ( LWORK.EQ.-1 )
*
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, 3*N-1 ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
*
      IF( INFO.EQ.0 ) THEN
         NB = ILAENV_( 1, 'DSYTRD', UPLO, N, -1, -1, -1 )
         LWKOPT = MAX( 1, ( NB+2 )*N )
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYEV ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      IF( N.EQ.1 ) THEN
         W( 1 ) = A( 1, 1 )
         WORK( 1 ) = 3
         IF( WANTZ )
     $      A( 1, 1 ) = ONE
         RETURN
      END IF
*
*     Get machine constants.
*
      SAFMIN = DLAMCH_( 'Safe minimum' )
      EPS = DLAMCH_( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = SQRT( BIGNUM )
*
*     Scale matrix to allowable range, if necessary.
*
      ANRM = DLANSY_( 'M', UPLO, N, A, LDA, WORK )
      ISCALE = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 )
     $   call dlascl_( UPLO, 0, 0, ONE, SIGMA, N, N, A, LDA, INFO )
*
*     Call DSYTRD to reduce symmetric matrix to tridiagonal form.
*
      INDE = 1
      INDTAU = INDE + N
      INDWRK = INDTAU + N
      LLWORK = LWORK - INDWRK + 1
      call dsytrd_( UPLO, N, A, LDA, W, WORK( INDE ), WORK( INDTAU ),
     $             WORK( INDWRK ), LLWORK, IINFO )
      LOPT = 2*N + Int(WORK( INDWRK ))
*
*     For eigenvalues only, call DSTERF.  For eigenvectors, first call
*     DORGTR to generate the orthogonal matrix, then call DSTEQR.
*
      IF( .NOT.WANTZ ) THEN
         call dsterf_( N, W, WORK( INDE ), INFO )
      ELSE
         call dorgtr_( UPLO, N, A, LDA, WORK( INDTAU ), WORK( INDWRK ),
     $                LLWORK, IINFO )
         call dsteqr_( JOBZ, N, W, WORK( INDE ), A, LDA, WORK( INDTAU ),
     $                INFO )
      END IF
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = N
         ELSE
            IMAX = INFO - 1
         END IF
         CALL DSCAL_( IMAX, ONE / SIGMA, W, 1 )
      END IF
*
*     Set WORK(1) to optimal workspace size.
*
      WORK( 1 ) = LWKOPT
*
      RETURN
*
*     End of DSYEV
*
      END
