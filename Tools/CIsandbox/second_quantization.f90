!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2013, Steven Vancoillie                                *
!***********************************************************************
MODULE SECOND_QUANTIZATION
  ! Implements determinants of k electrons in n (spin)orbitals as combinations
  ! (k,n) represented as bitstrings, using 32-bit integers.  Combinations are
  ! represented in lexicographic ordering, and their rank can be computed using
  ! a binomial table. The 32st bit is used for to indicate sign. All-one-bits
  ! are used to indicate annihilation (destroyed state). Excitation operators
  ! can be applied to the determinants and return the resulting determinant.

  ! Additionally, a wavefunction is represtented by a structure (derived type)
  ! consisting of determinant coefficients resulting from the direct product of
  ! separate alpha and beta determinants.

  ! Steven Vancoillie, November 2013, Lund

  USE ISO_FORTRAN_ENV, ONLY: INT64
  IMPLICIT NONE
  INTEGER, SAVE :: ONEBITS(0:255), RANKTBL(0:255,64)

CONTAINS

  RECURSIVE INTEGER FUNCTION GCD(I,J) RESULT (K)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: I, J
    IF (J.EQ.0) THEN
      K = I
    ELSE
      K = GCD(J,MOD(I,J))
    END IF
  END FUNCTION GCD

  INTEGER FUNCTION BINOM_COEF(K,N)
    ! COMPUTE BINOMIAL COEFFICIENT
    ! RETURNS #K-SUBSETS OUT OF N CHOICES
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: K, N
    INTEGER :: I, FRAC(2), DIV
    IF (K .GT. N) THEN
      BINOM_COEF = 0_INT64
    ELSE
      FRAC = 1_INT64
      DO I=1,K
        FRAC(1) = FRAC(1) * (N-K+I)
        FRAC(2) = FRAC(2) * I
        DIV = GCD(FRAC(1),FRAC(2))
        IF (DIV.GT.1) THEN
          FRAC = FRAC / DIV
        END IF
      END DO
      BINOM_COEF = FRAC(1)/FRAC(2)
    END IF
  END FUNCTION BINOM_COEF

  INTEGER FUNCTION LEX_INIT(K,N) RESULT(C)
    ! INITIALIZE THE FIRST COMBINATION
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: K, N
    IF (K .GT. N) THEN
      C = 0
    ELSE
      C = 2_INT64**K - 1
    END IF
  END FUNCTION LEX_INIT

  INTEGER FUNCTION LEX_NEXT(C) RESULT(D)
    ! GENERATE THE NEXT LEXICOGRAPHIC BIT COMBINATION
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: C
    INTEGER :: T, S
    T = IOR(C,C-1)+1
    S = IAND(NOT(T-1),T)-1
    D = IOR(T,ISHFT(S,-(TRAILZ(C)+1)))
  END FUNCTION LEX_NEXT

  SUBROUTINE RANK_INIT
    ! INITIALIZES A RANK TABLE INDEXED AS:
    ! 1 2   3    4        8
    !         1 1-8 1-16 1-24 ... 1-56
    !  0
    !  1
    !  ...
    !  255
    ! THE ROW INDEX IS A HASH TABLE USING THE
    ! BIT REPRESENTATION OF A BYTE SECTION.
    ! THE COLUMNS REFER TO THE POSITION OF
    ! A BYTE SECTION AND THE NUMBER OF POSSIBLE
    ! PRECEDING ONE-BITS
    INTEGER :: IROW, ICOL, IBYTE, IOFFSET
    INTEGER :: IBIT, IPOS
    INTEGER :: RANK

    DO IROW = 0, 255
      ONEBITS(IROW) = 0
      DO IBIT = 0, 7
        IF (BTEST(IROW,IBIT)) THEN
          ONEBITS(IROW) = ONEBITS(IROW) + 1
        END IF
      END DO
    END DO

    DO IROW = 0, 255
      ! COMPUTE RANK FOR FIRST BYTE
      RANK = 0
      IPOS = 0
      DO IBIT = 0, 7
        IF (BTEST(IROW,IBIT)) THEN
          IPOS = IPOS + 1
          RANK = RANK + BINOM_COEF(IPOS,IBIT)
        END IF
      END DO
      ICOL = 1
      RANKTBL(IROW,ICOL) = RANK
      ! COMPUTE RANK FOR SUBSEQUENT BYTES
      ! DEPENDING ON THE NUMBER OF PRECEDING
      ! ONE-BITS IN THE PREVIOUS BYTES
      DO IBYTE = 2, 4
        DO IOFFSET = 0, 8*(IBYTE-1)
          ICOL = ICOL + 1
          ! COMPUTE RANK
          RANK = 0
          IPOS = IOFFSET
          DO IBIT = 0, 7
            IF (BTEST(IROW,IBIT)) THEN
              IPOS = IPOS + 1
              RANK = RANK + BINOM_COEF(IPOS,IBIT+8*(IBYTE-1))
            END IF
          END DO
          RANKTBL(IROW,ICOL) = RANK
        END DO
      END DO
    END DO
  END SUBROUTINE RANK_INIT

  INTEGER FUNCTION RANK_(C)
    INTEGER, INTENT(IN) :: C
    INTEGER :: BYTE(4), ONES(4)
    RANK_ = 0
    IF (C.EQ.-1) RETURN
    BYTE(1) = IBITS(C, 0,8)
    BYTE(2) = IBITS(C, 8,8)
    BYTE(3) = IBITS(C,16,8)
    BYTE(4) = IBITS(C,24,6)
    ONES(1) = ONEBITS(BYTE(1))
    ONES(2) = ONEBITS(BYTE(2)) + ONES(1)
    ONES(3) = ONEBITS(BYTE(3)) + ONES(2)
    RANK_ = 1 + RANKTBL(BYTE(1), 1)     &
          & + RANKTBL(BYTE(2), 2+ONES(1)) &
          & + RANKTBL(BYTE(3),11+ONES(2)) &
          & + RANKTBL(BYTE(4),28+ONES(3))
  END FUNCTION RANK_

  INTEGER FUNCTION EX1(P,Q,C) RESULT(D)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: P, Q, C
    INTEGER :: T
    D = C
    IF (.NOT.BTEST(D,Q-1)) THEN
      D = -1
      RETURN
    END IF
    D = IBCLR(D,Q-1)
    IF (BTEST(D,P-1)) THEN
      D = -1
      RETURN
    END IF
    D = IBSET(D,P-1)
    IF(P.GT.Q) THEN
      T=IBITS(D,Q,P-Q-1)
    ELSE IF (P.LT.Q) THEN
      T=IBITS(D,P,Q-P-1)
    ELSE
      RETURN
    END IF
    DO WHILE (T.NE.0)
      D = IEOR(ISHFT(IAND(T,1),31),D)
      T = ISHFT(T,-1)
    END DO
  END FUNCTION EX1

  INTEGER FUNCTION ANN(P,C) RESULT(D)
    ! operates on determinant C with (â_p)
    ! and returns the resulting determinant with
    ! the sign in the highest bit (1 = negative)
    INTEGER, INTENT(IN) :: P, C
    INTEGER :: T
    IF (.NOT.BTEST(C,P-1)) THEN
      D = -1
      RETURN
    END IF
    D = IBCLR(C,P-1)
    T = IBITS(D,0,P-1)
    T = IEOR(T,ISHFT(T,-16))
    T = IEOR(T,ISHFT(T,-8))
    T = IEOR(T,ISHFT(T,-4))
    T = IAND(T,INT(B'1111'))
    T = IAND(ISHFT(INT(Z'6996'),-T),1)
    D = IEOR(ISHFT(T,31),D)
  END FUNCTION ANN

  INTEGER FUNCTION CRE(P,C) RESULT(D)
    ! operates on determinant C with (â†_p)
    ! and returns the resulting determinant with
    ! the sign in the highest bit (1 = negative)
    INTEGER, INTENT(IN) :: P, C
    INTEGER :: T
    IF (BTEST(C,P-1)) THEN
      D = -1
      RETURN
    END IF
    D = IBSET(C,P-1)
    T = IBITS(D,0,P-1)
    T = IEOR(T,ISHFT(T,-16))
    T = IEOR(T,ISHFT(T,-8))
    T = IEOR(T,ISHFT(T,-4))
    T = IAND(T,INT(B'1111'))
    T = IAND(ISHFT(INT(Z'6996'),-T),1)
    D = IEOR(ISHFT(T,31),D)
  END FUNCTION CRE

  INTEGER FUNCTION ANN2(P,Q,C) RESULT(D)
    ! operates on determinant C with (â_p â_q)
    ! and returns the resulting determinant with
    ! the sign in the highest bit (1 = negative)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: P, Q, C
    INTEGER :: T
    IF (.NOT.(BTEST(C,Q-1).AND.BTEST(C,P-1))) THEN
      D = -1
      RETURN
    END IF
    D = IBCLR(C,Q-1)
    D = IBCLR(D,P-1)
    IF(Q.LT.P) THEN
      T=IBITS(D,Q,P-Q-1)
    ELSE IF (P.LT.Q) THEN
      T=IBITS(D,P,Q-P-1)
      D=IEOR(ISHFT(1,31),D)
    ELSE
      D = -1
      RETURN
    END IF
    T = IEOR(T,ISHFT(T,-16))
    T = IEOR(T,ISHFT(T,-8))
    T = IEOR(T,ISHFT(T,-4))
    T = IAND(T,INT(B'1111'))
    T = IAND(ISHFT(INT(Z'6996'),-T),1)
    D = IEOR(ISHFT(T,31),D)
  END FUNCTION ANN2

  INTEGER FUNCTION CRE2(P,Q,C) RESULT(D)
    ! operates on determinant C with (â†_p â†_q)
    ! and returns the resulting determinant with
    ! the sign in the highest bit (1 = negative)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: P, Q, C
    INTEGER :: T
    IF ((BTEST(C,Q-1).OR.BTEST(C,P-1))) THEN
      D = -1
      RETURN
    END IF
    D = IBSET(C,Q-1)
    D = IBSET(D,P-1)
    IF(Q.LT.P) THEN
      T=IBITS(D,Q,P-Q-1)
      D=IEOR(ISHFT(1,31),D)
    ELSE IF (P.LT.Q) THEN
      T=IBITS(D,P,Q-P-1)
    ELSE
      D = -1
      RETURN
    END IF
    T = IEOR(T,ISHFT(T,-16))
    T = IEOR(T,ISHFT(T,-8))
    T = IEOR(T,ISHFT(T,-4))
    T = IAND(T,INT(B'1111'))
    T = IAND(ISHFT(INT(Z'6996'),-T),1)
    D = IEOR(ISHFT(T,31),D)
  END FUNCTION CRE2

  SUBROUTINE TEST_OPERATORS
    ! test routine for evaluation the creation and annihilation operators
    INTEGER :: NEL, NORB
    INTEGER :: NDET, IDET, DET, TMP
    INTEGER :: P

    ! test on 3in5
    NEL=3
    NORB=5

    NDET=BINOM_COEF(NEL,NORB)

    WRITE(*,*) 'Annihilation'
    DET=LEX_INIT(NEL,NORB)
    DO IDET=1,NDET
      WRITE(*,'(1X,B5.5)') IBITS(DET,0,5)
      DO P=1,NORB
        TMP=ANN(P,DET)
        IF (TMP.NE.-1) THEN
          WRITE (*,'(1X,2I3,1X,B5.5)') P, SIGN(1,TMP), IBITS(TMP,0,5)
        END IF
      END DO
      DET=LEX_NEXT(DET)
    END DO
    WRITE(*,*) 'Creation'
    DET=LEX_INIT(NEL,NORB)
    DO IDET=1,NDET
      WRITE(*,'(1X,B5.5)') IBITS(DET,0,5)
      DO P=1,NORB
        TMP=CRE(P,DET)
        IF (TMP.NE.-1) THEN
          WRITE (*,'(1X,2I3,1X,B5.5)') P, SIGN(1,TMP), IBITS(TMP,0,5)
        END IF
      END DO
      DET=LEX_NEXT(DET)
    END DO

  END SUBROUTINE TEST_OPERATORS

END MODULE SECOND_QUANTIZATION
