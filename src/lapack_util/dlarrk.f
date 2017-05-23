*> \brief \b DLARRK computes one eigenvalue of a symmetric tridiagonal matrix T to suitable accuracy.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DLARRK + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarrk.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarrk.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarrk.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLARRK( N, IW, GL, GU,
*                           D, E2, PIVMIN, RELTOL, W, WERR, INFO)
*
*       .. Scalar Arguments ..
*       INTEGER   INFO, IW, N
*       REAL*8              PIVMIN, RELTOL, GL, GU, W, WERR
*       ..
*       .. Array Arguments ..
*       REAL*8             D( * ), E2( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLARRK computes one eigenvalue of a symmetric tridiagonal
*> matrix T to suitable accuracy. This is an auxiliary code to be
*> called from DSTEMR.
*>
*> To avoid overflow, the matrix must be scaled so that its
*> largest element is no greater than overflow**(1/2) * underflow**(1/4) in absolute value, and for greatest
*> accuracy, it should not be much smaller than that.
*>
*> See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
*> Matrix", Report CS41, Computer Science Dept., Stanford
*> University, July 21, 1966.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the tridiagonal matrix T.  N >= 0.
*> \endverbatim
*>
*> \param[in] IW
*> \verbatim
*>          IW is INTEGER
*>          The index of the eigenvalues to be returned.
*> \endverbatim
*>
*> \param[in] GL
*> \verbatim
*>          GL is REAL*8
*> \endverbatim
*>
*> \param[in] GU
*> \verbatim
*>          GU is REAL*8
*>          An upper and a lower bound on the eigenvalue.
*> \endverbatim
*>
*> \param[in] D
*> \verbatim
*>          D is REAL*8           array, dimension (N)
*>          The n diagonal elements of the tridiagonal matrix T.
*> \endverbatim
*>
*> \param[in] E2
*> \verbatim
*>          E2 is REAL*8           array, dimension (N-1)
*>          The (n-1) squared off-diagonal elements of the tridiagonal matrix T.
*> \endverbatim
*>
*> \param[in] PIVMIN
*> \verbatim
*>          PIVMIN is REAL*8
*>          The minimum pivot allowed in the Sturm sequence for T.
*> \endverbatim
*>
*> \param[in] RELTOL
*> \verbatim
*>          RELTOL is REAL*8
*>          The minimum relative width of an interval.  When an interval
*>          is narrower than RELTOL times the larger (in
*>          magnitude) endpoint, then it is considered to be
*>          sufficiently small, i.e., converged.  Note: this should
*>          always be at least radix*machine epsilon.
*> \endverbatim
*>
*> \param[out] W
*> \verbatim
*>          W is REAL*8
*> \endverbatim
*>
*> \param[out] WERR
*> \verbatim
*>          WERR is REAL*8
*>          The error bound on the corresponding eigenvalue approximation
*>          in W.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:       Eigenvalue converged
*>          = -1:      Eigenvalue did NOT converge
*> \endverbatim
*
*> \par Internal Parameters:
*  =========================
*>
*> \verbatim
*>  FUDGE   REAL*8          , default = 2
*>          A "fudge factor" to widen the Gershgorin intervals.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date September 2012
*
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      SUBROUTINE DLARRK( N, IW, GL, GU,
     $                    D, E2, PIVMIN, RELTOL, W, WERR, INFO)
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      INTEGER   INFO, IW, N
      REAL*8              PIVMIN, RELTOL, GL, GU, W, WERR
*     ..
*     .. Array Arguments ..
      REAL*8             D( * ), E2( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL*8             FUDGE, HALF, TWO, ZERO
      PARAMETER          ( HALF = 0.5D0, TWO = 2.0D0,
     $                     FUDGE = TWO, ZERO = 0.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER   I, IT, ITMAX, NEGCNT
      REAL*8             ATOLI, EPS, LEFT, MID, RIGHT, RTOLI, TMP1,
     $                   TMP2, TNORM
*     ..
*     .. External Functions ..
      REAL*8             DLAMCH
      EXTERNAL   DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX
*     ..
*     .. Executable Statements ..
*
*     Get machine constants
      EPS = DLAMCH( 'P' )

      TNORM = MAX( ABS( GL ), ABS( GU ) )
      RTOLI = RELTOL
      ATOLI = FUDGE*TWO*PIVMIN

      ITMAX = INT( ( LOG( TNORM+PIVMIN )-LOG( PIVMIN ) ) /
     $           LOG( TWO ) ) + 2

      INFO = -1

      LEFT = GL - FUDGE*TNORM*EPS*N - FUDGE*TWO*PIVMIN
      RIGHT = GU + FUDGE*TNORM*EPS*N + FUDGE*TWO*PIVMIN
      IT = 0

 10   CONTINUE
*
*     Check if interval converged or maximum number of iterations reached
*
      TMP1 = ABS( RIGHT - LEFT )
      TMP2 = MAX( ABS(RIGHT), ABS(LEFT) )
      IF( TMP1.LT.MAX( ATOLI, PIVMIN, RTOLI*TMP2 ) ) THEN
         INFO = 0
         GOTO 30
      ENDIF
      IF(IT.GT.ITMAX)
     $   GOTO 30

*
*     Count number of negative pivots for mid-point
*
      IT = IT + 1
      MID = HALF * (LEFT + RIGHT)
      NEGCNT = 0
      TMP1 = D( 1 ) - MID
      IF( ABS( TMP1 ).LT.PIVMIN )
     $   TMP1 = -PIVMIN
      IF( TMP1.LE.ZERO )
     $   NEGCNT = NEGCNT + 1
*
      DO 20 I = 2, N
         TMP1 = D( I ) - E2( I-1 ) / TMP1 - MID
         IF( ABS( TMP1 ).LT.PIVMIN )
     $      TMP1 = -PIVMIN
         IF( TMP1.LE.ZERO )
     $      NEGCNT = NEGCNT + 1
 20   CONTINUE

      IF(NEGCNT.GE.IW) THEN
         RIGHT = MID
      ELSE
         LEFT = MID
      ENDIF
      GOTO 10

 30   CONTINUE
*
*     Converged or maximum number of iterations reached
*
      W = HALF * (LEFT + RIGHT)
      WERR = HALF * ABS( RIGHT - LEFT )

      RETURN
*
*     End of DLARRK
*
      END SUBROUTINE
