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
C     Expand IMS component of ICS(NLEV) in determinants using the
C     procedure from Shavitt in "The Unitary Group", "Lecture Notes in
C     Chemistry" Vol. 22, pp. 55.
      SUBROUTINE EXPCSF (ICS, NLEV, IMS, LEX)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ICS(NLEV), LEX(NLEV)
      CHARACTER(256) LINE
      CHARACTER(6) STRING
      DIMENSION ICOEF(2)
      LOGICAL LPHASE, LAST
C     find number of singly occupied orbitals
      NSOMO = 0
      DO I=1,NLEV
        IF (ICS(I).EQ.1.OR.ICS(I).EQ.2) NSOMO=NSOMO+1
      ENDDO
      NALPHA=(NSOMO+IMS)/2
      CALL INIT_LEX (NSOMO, NALPHA, LEX)
C     Loop over possible determinants
      LAST = .FALSE.
      DO WHILE (.NOT.LAST)
        ICOEF(1)=1
        ICOEF(2)=1
        LPHASE=.TRUE.
        IALPHA=0
        IBETA=0
        ISOMO=0
        ILEX=1
        IA=0
        IB=0
        LINE=' '
        K=26
        WRITE(LINE(K:K),'(A)') '|'
        K=K+NLEV+1
        WRITE(LINE(K:K),'(A)') '|'
        K=26
        DO I=1,NLEV
          K=K+1
          IF (ICS(I).EQ.0) THEN
            LINE(K:K)='0'
          ELSE IF (ICS(I).EQ.1) THEN
            IB=IB+1
            ISOMO=ISOMO+1
            IF (ISOMO.EQ.LEX(ILEX)) THEN
              ILEX=ILEX+1
              ICOEF(1)=ICOEF(1)*(IA+IB-IBETA)
              IALPHA=IALPHA+1
              LINE(K:K)='a'
            ELSE
              ICOEF(1)=ICOEF(1)*(IA+IB-IALPHA)
              IBETA=IBETA+1
              LINE(K:K)='b'
            ENDIF
            ICOEF(2)=ICOEF(2)*IB
          ELSE IF (ICS(I).EQ.2) THEN
            IA=IA+1
            IB=IB-1
            ISOMO=ISOMO+1
            IF (ISOMO.EQ.LEX(ILEX)) THEN
              ILEX=ILEX+1
              ICOEF(1)=ICOEF(1)*(IBETA-IA+1)
              IF(MOD(IB,2).EQ.0) LPHASE=.NOT.LPHASE
              IALPHA=IALPHA+1
              LINE(K:K)='a'
            ELSE
              ICOEF(1)=ICOEF(1)*(IALPHA-IA+1)
              IF(MOD(IB,2).NE.0) LPHASE=.NOT.LPHASE
              IBETA=IBETA+1
              LINE(K:K)='b'
            ENDIF
            ICOEF(2)=ICOEF(2)*(IB+2)
          ELSE
            IA=IA+1
            IALPHA=IALPHA+1
            IBETA=IBETA+1
            IF (MOD(IB,2).NE.0) LPHASE=.NOT.LPHASE
            LINE(K:K)='2'
          ENDIF
          IF (ICOEF(1).EQ.0) EXIT
        ENDDO
        IF (ICOEF(1).NE.0) THEN
C     If non-zero coefficient, print the thing
          CALL SIMPLIFY(ICOEF)
          IF (LPHASE) THEN
            WRITE (LINE(1:9),'(2X,A7)') '+ sqrt('
          ELSE
            WRITE (LINE(1:9),'(2X,A7)') '- sqrt('
          ENDIF
          WRITE(LINE(10:16),'(I6,A1)') ICOEF(1), '/'
          WRITE(STRING(1:6),'(I6)') ICOEF(2)
          J=17
          DO I=1,6
            IF (STRING(I:I).NE.' ') THEN
              LINE(J:J)=STRING(I:I)
              J=J+1
            ENDIF
          ENDDO
          WRITE(LINE(23:23),'(A1)') ')'
          WRITE(6,*) LINE(1:K+1)
C     End of print-out
        ENDIF
C     Get the next determinant
        CALL LEX_ITER (NSOMO, NALPHA, LEX, LAST)
      ENDDO
      END

      SUBROUTINE INIT_LEX (N, K, LEX)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION LEX(K)
      DO I=1,K
        LEX(I)=I
      ENDDO
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(N)
      END

      SUBROUTINE LEX_ITER (N, K, LEX, LAST)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION LEX(K)
      LOGICAL LAST
      I = K
C     Get the first position to be updated
      DO WHILE ((I.GT.0).AND.(LEX(I).EQ.N-K+I))
        I=I-1
      ENDDO
C     If still remaining combinations, update and
C     reset all higher positions to lexicographic order
      IF (I.GT.0) THEN
        LEX(I)=LEX(I)+1
        DO J=1,K-I
          LEX(I+J)=LEX(I)+J
        ENDDO
C     Else, quit finding combinations
      ELSE
        LAST = .TRUE.
      ENDIF
      END

      SUBROUTINE SIMPLIFY (FRAC)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER A, B, FRAC
      DIMENSION FRAC(2)
      IF (FRAC(1).EQ.0) RETURN
C     Find GCD of numerator and denominator
      A=FRAC(1)
      B=FRAC(2)
 60   IF (B.NE.0) THEN
        ITEMP=B
        B=MOD(A,B)
        A=ITEMP
        GOTO 60
      ENDIF
      FRAC(1)=FRAC(1)/A
      FRAC(2)=FRAC(2)/A
      END
