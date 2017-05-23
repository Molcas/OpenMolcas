      SUBROUTINE DLAMC1( BETA, T, RND, IEEE1 )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
*     ..
*
*  Purpose
*  =======
*
*  DLAMC1 determines the machine parameters given by BETA, T, RND, and
*  IEEE1.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  IEEE1   (output) LOGICAL
*          Specifies whether rounding appears to be done in the IEEE
*          'round to nearest' style.
*
*  Further Details
*  ===============
*
*  The routine is based on the routine  ENVRON  by Malcolm and
*  incorporates suggestions by Gentleman and Marovich. See
*
*     Malcolm M. A. (1972) Algorithms to reveal properties of
*        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
*
*     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
*        that reveal properties of floating point arithmetic units.
*        Comms. of the ACM, 17, 276-277.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      REAL*8             A, B, C, F, HALF, ONE, QTR, SAVEC, T1, T2
*     ..
*     .. External Functions ..
      REAL*8             DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Save statement ..
      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ONE = 1
         HALF = ONE/2
*
*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
*        IEEE1, T and RND.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are  stored and not held in registers,  or
*        are not affected by optimizers.
*
*        Compute  a = 2.0**m  with the  smallest positive integer m such
*        that
*
*           fl( a + 1.0 ) = a.
*
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   10    CONTINUE
         IF( C.EQ.ONE ) THEN
            A = 2*A
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 10
         END IF
*+       END WHILE
*
*        Now compute  b = 2.0**m  with the smallest positive integer m
*        such that
*
*           fl( a + b ) .gt. a.
*
         B = 1
         C = DLAMC3( A, B )
         C = DLAMC3( C, -A )
*
*+       WHILE( C.EQ.A )LOOP
   20    CONTINUE
         IF( C.LE.HALF ) THEN
            B = 2*B
            C = DLAMC3( A, B )
            C = DLAMC3( C, -A )
            GO TO 20
         END IF
*+       END WHILE
         C = DLAMC3( A, B )
*
*        Now compute the base.  a and c  are neighbouring floating point
*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
*        their difference is beta. Adding 0.25 to c is to ensure that it
*        is truncated to beta and not ( beta - 1 ).
*
         QTR = ONE / 4
         SAVEC = C
         C = DLAMC3( C, -A )
         LBETA = Int( C + QTR )
*
*        Now determine whether rounding or chopping occurs,  by adding a
*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
*
         B = LBETA
         F = DLAMC3( B / 2, -B / 100 )
         C = DLAMC3( F, A )
         IF( C.EQ.A ) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = DLAMC3( B / 2, B / 100 )
         C = DLAMC3( F, A )
         IF( ( LRND ) .AND. ( C.EQ.A ) )
     $      LRND = .FALSE.
*
*        Try and decide whether rounding is done in the  IEEE  'round to
*        nearest' style. B/2 is half a unit in the last place of the two
*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
*        A, but adding B/2 to SAVEC should change SAVEC.
*
         T1 = DLAMC3( B / 2, A )
         T2 = DLAMC3( B / 2, SAVEC )
         LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND
*
*        Now find  the  mantissa, t.  It should  be the  integer part of
*        log to the base beta of a,  however it is safer to determine  t
*        by powering.  So we find t as the smallest positive integer for
*        which
*
*           fl( beta**t + 1.0 ) = 1.0.
*
         LT = 0
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   30    CONTINUE
         IF( C.EQ.ONE ) THEN
            LT = LT + 1
            A = A*LBETA
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 30
         END IF
*+       END WHILE
*
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      RETURN
*
*     End of DLAMC1
*
      END SUBROUTINE
corg      SUBROUTINE DLAMC1( BETA, T, RND, IEEE1 )
corg*
corg*  -- LAPACK auxiliary routine (version 3.0) --
corg*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
corg*     Courant Institute, Argonne National Lab, and Rice University
corg*     October 31, 1992
corg*
corg*     .. Scalar Arguments ..
corg      LOGICAL            IEEE1, RND
corg      INTEGER            BETA, T
corg*     ..
corg*
corg*  Purpose
corg*  =======
corg*
corg*  DLAMC1 determines the machine parameters given by BETA, T, RND, and
corg*  IEEE1.
corg*
corg*  Arguments
corg*  =========
corg*
corg*  BETA    (output) INTEGER
corg*          The base of the machine.
corg*
corg*  T       (output) INTEGER
corg*          The number of ( BETA ) digits in the mantissa.
corg*
corg*  RND     (output) LOGICAL
corg*          Specifies whether proper rounding  ( RND = .TRUE. )  or
corg*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
corg*          be a reliable guide to the way in which the machine performs
corg*          its arithmetic.
corg*
corg*  IEEE1   (output) LOGICAL
corg*          Specifies whether rounding appears to be done in the IEEE
corg*          'round to nearest' style.
corg*
corg*  Further Details
corg*  ===============
corg*
corg*  The routine is based on the routine  ENVRON  by Malcolm and
corg*  incorporates suggestions by Gentleman and Marovich. See
corg*
corg*     Malcolm M. A. (1972) Algorithms to reveal properties of
corg*        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
corg*
corg*     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
corg*        that reveal properties of floating point arithmetic units.
corg*        Comms. of the ACM, 17, 276-277.
corg*
corg* =====================================================================
corg*
corg*     .. Local Scalars ..
corg      LOGICAL            FIRST, LIEEE1, LRND
corg      INTEGER            LBETA, LT
corg      REAL*8             A, B, C, F, ONE, QTR, SAVEC, T1, T2
corg*     ..
corg*     .. External Functions ..
corg      REAL*8             DLAMC3
corg      EXTERNAL           DLAMC3
corg*     ..
corg*     .. Save statement ..
corg      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
corg*     ..
corg*     .. Data statements ..
corg      DATA               FIRST / .TRUE. /
corg*     ..
corg*     .. Executable Statements ..
corg*
corg      IF( FIRST ) THEN
corg         FIRST = .FALSE.
corg         ONE = 1
corg*
corg*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
corg*        IEEE1, T and RND.
corg*
corg*        Throughout this routine  we use the function  DLAMC3  to ensure
corg*        that relevant values are  stored and not held in registers,  or
corg*        are not affected by optimizers.
corg*
corg*        Compute  a = 2.0**m  with the  smallest positive integer m such
corg*        that
corg*
corg*           fl( a + 1.0 ) = a.
corg*
corg         A = 1
corg         C = 1
corg*
corg*+       WHILE( C.EQ.ONE )LOOP
corg   10    CONTINUE
corg         IF( C.EQ.ONE ) THEN
corg            A = 2*A
corg            C = DLAMC3( A, ONE )
corg            C = DLAMC3( C, -A )
corg            GO TO 10
corg         END IF
corg*+       END WHILE
corg*
corg*        Now compute  b = 2.0**m  with the smallest positive integer m
corg*        such that
corg*
corg*           fl( a + b ) .gt. a.
corg*
corg         B = 1
corg         C = DLAMC3( A, B )
corg*
corg*+       WHILE( C.EQ.A )LOOP
corg   20    CONTINUE
corg         IF( C.EQ.A ) THEN
corg            B = 2*B
corg            C = DLAMC3( A, B )
corg            GO TO 20
corg         END IF
corg*+       END WHILE
corg*
corg*        Now compute the base.  a and c  are neighbouring floating point
corg*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
corg*        their difference is beta. Adding 0.25 to c is to ensure that it
corg*        is truncated to beta and not ( beta - 1 ).
corg*
corg         QTR = ONE / 4
corg         SAVEC = C
corg         C = DLAMC3( C, -A )
corg         LBETA = C + QTR
corg*
corg*        Now determine whether rounding or chopping occurs,  by adding a
corg*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
corg*
corg         B = LBETA
corg         F = DLAMC3( B / 2, -B / 100 )
corg         C = DLAMC3( F, A )
corg         IF( C.EQ.A ) THEN
corg            LRND = .TRUE.
corg         ELSE
corg            LRND = .FALSE.
corg         END IF
corg         F = DLAMC3( B / 2, B / 100 )
corg         C = DLAMC3( F, A )
corg         IF( ( LRND ) .AND. ( C.EQ.A ) )
corg     $      LRND = .FALSE.
corg*
corg*        Try and decide whether rounding is done in the  IEEE  'round to
corg*        nearest' style. B/2 is half a unit in the last place of the two
corg*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
corg*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
corg*        A, but adding B/2 to SAVEC should change SAVEC.
corg*
corg         T1 = DLAMC3( B / 2, A )
corg         T2 = DLAMC3( B / 2, SAVEC )
corg         LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND
corg*
corg*        Now find  the  mantissa, t.  It should  be the  integer part of
corg*        log to the base beta of a,  however it is safer to determine  t
corg*        by powering.  So we find t as the smallest positive integer for
corg*        which
corg*
corg*           fl( beta**t + 1.0 ) = 1.0.
corg*
corg         LT = 0
corg         A = 1
corg         C = 1
corg*
corg*+       WHILE( C.EQ.ONE )LOOP
corg   30    CONTINUE
corg         IF( C.EQ.ONE ) THEN
corg            LT = LT + 1
corg            A = A*LBETA
corg            C = DLAMC3( A, ONE )
corg            C = DLAMC3( C, -A )
corg            GO TO 30
corg         END IF
corg*+       END WHILE
corg*
corg      END IF
corg*
corg      BETA = LBETA
corg      T = LT
corg      RND = LRND
corg      IEEE1 = LIEEE1
corg      RETURN
corg*
corg*     End of DLAMC1
corg*
corg      END
*
************************************************************************
*
