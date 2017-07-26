#ifndef _HAVE_EXTRA_

      SUBROUTINE EIGEN(A,R,N,MV,MFKR)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*),R(*)
      CALL Unused_Real_Array(A)
      CALL Unused_Real_Array(R)
      CALL Unused_Integer(N)
      CALL Unused_Integer(MV)
      CALL Unused_Integer(MFKR)
      END SUBROUTINE EIGEN

#endif
