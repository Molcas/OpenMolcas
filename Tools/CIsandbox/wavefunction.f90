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
MODULE WAVEFUNCTION
  ! Implements a wavefunction represtented by a structure (derived type)
  ! consisting of determinant coefficients resulting from the direct product of
  ! separate alpha and beta determinant strings.

  ! Steven Vancoillie, November 2013, Lund

  USE ISO_FORTRAN_ENV
  USE SECOND_QUANTIZATION
  IMPLICIT NONE

  ! determinant derived type
  TYPE DET
    INTEGER :: ALFA, BETA
  END type DET

  ! wavefunction derived type
  TYPE WFN
    INTEGER :: NEL, NORB, MULT
    INTEGER :: NELA, NDETA
    INTEGER :: NELB, NDETB
    INTEGER :: NDET
    REAL(REAL64), ALLOCATABLE :: COEF(:,:)
  END TYPE WFN

CONTAINS

  SUBROUTINE DET_PRINT(D,NORB)
    IMPLICIT NONE
    TYPE(DET), INTENT(IN) :: D
    INTEGER, INTENT(IN) :: NORB
    CHARACTER(LEN=64) FMT
    WRITE(FMT,'("(2(X,B",I0,".",I0"))")') NORB, NORB
    WRITE(*,FMT) D%ALFA, D%BETA
  END SUBROUTINE DET_PRINT

  SUBROUTINE WFN_INIT(PSI,NEL,NORB,MULT)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NEL, NORB, MULT
    TYPE(WFN) :: PSI
    ! set basic properties
    PSI%NEL = NEL
    PSI%NORB = NORB
    PSI%MULT = MULT
    ! compute alpha/beta subsets
    PSI%NELA = (NEL+(MULT-1))/2
    PSI%NELB = (NEL-(MULT-1))/2
    IF (PSI%NELA+PSI%NELB.NE.NEL) THEN
      STOP 'WFN_INIT: NELA+NELB != NEL'
    END IF
    PSI%NDETA = BINOM_COEF(PSI%NELA,PSI%NORB)
    PSI%NDETB = BINOM_COEF(PSI%NELB,PSI%NORB)
    PSI%NDET=PSI%NDETA*PSI%NDETB
    ALLOCATE(PSI%COEF(PSI%NDETA,PSI%NDETB))
  END SUBROUTINE WFN_INIT

  SUBROUTINE WFN_FREE(PSI)
    IMPLICIT NONE
    TYPE(WFN) :: PSI
    DEALLOCATE(PSI%COEF)
  END SUBROUTINE WFN_FREE

  SUBROUTINE WFN_NORMALIZE(PSI)
    IMPLICIT NONE
    TYPE(WFN), INTENT(INOUT) :: PSI
    INTEGER :: IDETA, IDETB
    REAL(REAL64) :: NF
    ! compute normalization factor
    NF=0.0D0
    DO IDETB=1,PSI%NDETB
      DO IDETA=1,PSI%NDETA
        NF=NF+PSI%COEF(IDETA,IDETB)**2
      END DO
    END DO
    NF=SQRT(NF)
    PSI%COEF=PSI%COEF/NF
  END SUBROUTINE WFN_NORMALIZE

  SUBROUTINE WFN_PRINT(PSI,THR)
    IMPLICIT NONE
    CHARACTER(LEN=256) FMT
    TYPE(WFN), INTENT(IN) :: PSI
    INTEGER :: IA, IB
    INTEGER :: DETA, DETB
    REAL(REAL64) :: THR
    INTEGER :: NBITS
    NBITS=MAX(PSI%NORB,5)
    WRITE(*,'(X,A)') 'WAVEFUNCTION SPECS:'
    WRITE(*,'(X,A)') '==================='
    WRITE(*,'(X,A,I12)') '#ELECTRONS    = ', PSI%NEL
    WRITE(*,'(X,A,I12)') '#ORBITALS     = ', PSI%NORB
    WRITE(*,'(X,A,I12)') '#MULTIPLICITY = ', PSI%MULT
    WRITE(*,'(X,A,I12)') '#DET (ALPHA)  = ', PSI%NDETA
    WRITE(*,'(X,A,I12)') '#DET (BETA)   = ', PSI%NDETB
    WRITE(*,'(X,A,I12)') '#DET (TOTAL)  = ', PSI%NDETA*PSI%NDETB
    WRITE(*,*)
    WRITE(FMT,'("(A18,X,2(2X,A",I0,"))")') NBITS
    WRITE(*,FMT) 'COEFFICIENTS','ALFA','BETA'
    WRITE(FMT,'("(F18.14,X,2(2X,B",I0,".",I0,"))")') NBITS, PSI%NORB
    DETB = LEX_INIT(PSI%NELB,PSI%NORB)
    DO IB = 1, PSI%NDETB
      DETA = LEX_INIT(PSI%NELA,PSI%NORB)
      DO IA = 1, PSI%NDETA
        IF (ABS(PSI%COEF(IA,IB)).GE.THR) WRITE(*,FMT) PSI%COEF(IA,IB), DETA, DETB
        DETA = LEX_NEXT(DETA)
      END DO
      DETB = LEX_NEXT(DETB)
    END DO
    WRITE(*,*)
  END SUBROUTINE WFN_PRINT

  SUBROUTINE DET_EX1(FACT,P,Q,DETA,DETB,TAU)
    ! computes |TAU> = |TAU> + FACT * Ê_pq |DETA,DETB>
    USE SECOND_QUANTIZATION
    IMPLICIT NONE
    
    REAL(REAL64),INTENT(IN) :: FACT
    INTEGER, INTENT(IN) :: P, Q
    INTEGER, INTENT(IN) :: DETA, DETB
    TYPE(WFN), INTENT(INOUT) :: TAU

    INTEGER :: IA, IB
    INTEGER :: TMPA, TMPB
    REAL(REAL64) :: FA, FB
    
    ! compute E_(pa)(qa)
    TMPA=EX1(P,Q,DETA)
    IF (TMPA.NE.-1) THEN
      IA=RANK_(TMPA)
      IB=RANK_(DETB)
      FA=SIGN(1,TMPA)
      TAU%COEF(IA,IB)=TAU%COEF(IA,IB)+FA*FACT
    END IF
    ! compute E_(pb)(qb)
    TMPB=EX1(P,Q,DETB)
    IF (TMPB.NE.-1) THEN
      IA=RANK_(DETA)
      IB=RANK_(TMPB)
      FB=SIGN(1,TMPB)
      TAU%COEF(IA,IB)=TAU%COEF(IA,IB)+FB*FACT
    END IF

  END SUBROUTINE DET_EX1

  SUBROUTINE WFN_EX1(FACT,P,Q,SGM,TAU)
    USE SECOND_QUANTIZATION
    IMPLICIT NONE
    REAL(REAL64),INTENT(IN) :: FACT
    INTEGER, INTENT(IN) :: P, Q
    TYPE(WFN), INTENT(IN) :: SGM
    TYPE(WFN), INTENT(INOUT) :: TAU

    INTEGER :: IA, IB, DETA, DETB
    REAL(REAL64) :: CPQ

    DETB=LEX_INIT(SGM%NELB,SGM%NORB)
    DO IB=1,SGM%NDETB
      DETA=LEX_INIT(SGM%NELA,SGM%NORB)
      DO IA=1,SGM%NDETA

        CPQ=FACT*SGM%COEF(IA,IB)
        CALL DET_EX1(CPQ,P,Q,DETA,DETB,TAU)

        DETA=LEX_NEXT(DETA)
      END DO
      DETB=LEX_NEXT(DETB)
    END DO

  END SUBROUTINE WFN_EX1

  SUBROUTINE DET_EX2(FACT,P,Q,R,S,DETA,DETB,TAU)
    ! computes |TAU> = |TAU> + FACT * Ê_pqrs |DETA,DETB>
    USE SECOND_QUANTIZATION
    IMPLICIT NONE

    REAL(REAL64),INTENT(IN) :: FACT
    INTEGER, INTENT(IN) :: P, Q, R, S
    INTEGER, INTENT(IN) :: DETA, DETB
    TYPE(WFN), INTENT(INOUT) :: TAU

    INTEGER :: IA, IB
    INTEGER :: TMPA, TMPB
    REAL(REAL64) :: FA, FB

    ! compute E_(pa)(qa)(ra)(sa)
    TMPA=ANN2(S,Q,DETA)
    IF (TMPA.NE.-1) THEN
      TMPA=CRE2(P,R,TMPA)
      IF (TMPA.NE.-1) THEN
        IA=RANK_(TMPA)
        IB=RANK_(DETB)
        FA=SIGN(1,TMPA)
        TAU%COEF(IA,IB)=TAU%COEF(IA,IB)+FA*FACT
      END IF
    END IF
    ! compute E_(pb)(qb)(rb)(sb)
    TMPB=ANN2(S,Q,DETB)
    IF (TMPB.NE.-1) THEN
      TMPB=CRE2(P,R,TMPB)
      IF (TMPB.NE.-1) THEN
        IA=RANK_(DETA)
        IB=RANK_(TMPB)
        FB=SIGN(1,TMPB)
        TAU%COEF(IA,IB)=TAU%COEF(IA,IB)+FB*FACT
      END IF
    END IF
    ! compute E_(pa)(qa)(rb)(sb)
    TMPA=ANN(Q,DETA)
    TMPB=ANN(S,DETB)
    IF (TMPA.NE.-1.AND.TMPB.NE.-1) THEN
      TMPA=CRE(P,TMPA)
      TMPB=CRE(R,TMPB)
      IF (TMPA.NE.-1.AND.TMPB.NE.-1) THEN
        IA=RANK_(TMPA)
        IB=RANK_(TMPB)
        FA=SIGN(1,TMPA)
        FB=SIGN(1,TMPB)
        TAU%COEF(IA,IB)=TAU%COEF(IA,IB)+FA*FB*FACT
      END IF
    END IF
    ! compute E_(pb)(qb)(ra)(sa)
    TMPA=ANN(S,DETA)
    TMPB=ANN(Q,DETB)
    IF (TMPA.NE.-1.AND.TMPB.NE.-1) THEN
      TMPA=CRE(R,TMPA)
      TMPB=CRE(P,TMPB)
      IF (TMPA.NE.-1.AND.TMPB.NE.-1) THEN
        IA=RANK_(TMPA)
        IB=RANK_(TMPB)
        FA=SIGN(1,TMPA)
        FB=SIGN(1,TMPB)
        TAU%COEF(IA,IB)=TAU%COEF(IA,IB)+FA*FB*FACT
      END IF
    END IF

  END SUBROUTINE DET_EX2

  SUBROUTINE WFN_COPY0(NI,SGM,TAU)
    ! sets |TAU> = |SGM>, where SGM excludes the inactives
    USE SECOND_QUANTIZATION
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NI
    TYPE(WFN), INTENT(IN) :: SGM
    TYPE(WFN), INTENT(INOUT) :: TAU
    INTEGER :: DETA, DETB, IA, IB
    INTEGER :: TAUA, TAUB
    INTEGER :: INACTIVE

    TAU%COEF=0.0D0

    INACTIVE=2**NI-1
    DETB=LEX_INIT(SGM%NELB,SGM%NORB)
    DO IB=1,SGM%NDETB
      TAUB=ISHFT(DETB,NI)+INACTIVE
      DETA=LEX_INIT(SGM%NELA,SGM%NORB)
      DO IA=1,SGM%NDETA
        TAUA=ISHFT(DETA,NI)+INACTIVE
        TAU%COEF(RANK_(TAUA),RANK_(TAUB))=SGM%COEF(IA,IB)
        DETA=LEX_NEXT(DETA)
      END DO
      DETB=LEX_NEXT(DETB)
    END DO
  END SUBROUTINE WFN_COPY0

  SUBROUTINE WFN_PK(NI,NA,NS,SGM,TAU)
    ! computes |TAU> = |TAU> + P_K |SGM>
    USE SECOND_QUANTIZATION
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NI, NA, NS
    TYPE(WFN), INTENT(IN) :: SGM
    TYPE(WFN), INTENT(INOUT) :: TAU
    INTEGER :: IA, IB, DETA, DETB, MASK

    ! keep only determinants from |SGM> which have all bits set in the lower NI
    ! part and no bits set in the upper NS part.
    MASK=2**NI-1+ISHFT(2**NS-1,NI+NA)
    DETB=LEX_INIT(SGM%NELB,SGM%NORB)
    DO IB=1,SGM%NDETB
      IF (IAND(DETB,MASK).EQ.2**NI-1) THEN
        DETA=LEX_INIT(SGM%NELA,SGM%NORB)
        DO IA=1,SGM%NDETA
          IF (IAND(DETA,MASK).EQ.2**NI-1) THEN
            TAU%COEF(IA,IB)=TAU%COEF(IA,IB)+SGM%COEF(IA,IB)
          END IF
          DETA=LEX_NEXT(DETA)
        END DO
      END IF
      DETB=LEX_NEXT(DETB)
    END DO
  END SUBROUTINE WFN_PK
  
  SUBROUTINE WFN_PSD(NI,NA,NS,SGM,TAU)
    ! computes |TAU> = |TAU> + P_SD |SGM>
    USE SECOND_QUANTIZATION
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NI, NA, NS
    TYPE(WFN), INTENT(IN) :: SGM
    TYPE(WFN), INTENT(INOUT) :: TAU
    INTEGER :: IA, IB, DETA, DETB, IMASK, SMASK
    INTEGER :: NIA, NSA, NIB, NSB
    INTEGER :: ICOUNT

    ! keep only determinants from |SGM> which have at least 2*NI-2 bits set in
    ! the lower NI part and at most two bits set in the upper NS part, and
    ! remove any non-excited determinants.

    ICOUNT=0
    IMASK=2**NI-1
    SMASK=ISHFT(2**NS-1,NI+NA)
    DETB=LEX_INIT(SGM%NELB,SGM%NORB)
    DO IB=1,SGM%NDETB
      NIB=ONEBITS(IAND(DETB,IMASK))
      NSB=ONEBITS(IAND(DETB,SMASK))
      DETA=LEX_INIT(SGM%NELA,SGM%NORB)
      DO IA=1,SGM%NDETA
        NIA=ONEBITS(IAND(DETA,IMASK))
        NSA=ONEBITS(IAND(DETA,SMASK))
        IF (NIB+NIA.GE.2*NI-2.AND.NSB+NSA.LE.2) THEN
          IF(.NOT.(NIB+NIA.EQ.2*NI.AND.NSB+NSA.EQ.0)) THEN
            TAU%COEF(IA,IB)=TAU%COEF(IA,IB)+SGM%COEF(IA,IB)
            ICOUNT=ICOUNT+1
          END IF
        END IF
        DETA=LEX_NEXT(DETA)
      END DO
      DETB=LEX_NEXT(DETB)
    END DO
  END SUBROUTINE WFN_PSD
  
  SUBROUTINE WFN_ORBROT(PSI,OVEC)
    ! Transforms the CI coefficients to match the orbital rotation given by the
    ! transformation matrix U. The dimensions of U should match those of the
    ! number of orbitals in the wavefunction (of course).
    !
    ! This one-orbital-at-a-time algorithm was adopted from the TRACI_RPT2
    ! subroutine in Molcas, and more info can be found in part IV (on page 4505)
    ! of the publication in Phys. Rev. E, vol. 52 (1995), p. 4499-4508.

    IMPLICIT NONE

    TYPE(WFN), INTENT(INOUT) :: PSI
    REAL(REAL64), INTENT(IN) :: OVEC(:,:)

    TYPE(WFN) :: TAU
    REAL(REAL64), ALLOCATABLE :: U(:,:), T(:)

    INTEGER :: NA, I, J, K
    REAL(REAL64) :: FACT

    NA=SIZE(OVEC,1)
    IF (PSI%NORB.NE.NA) STOP 'WFN_ORBROT: NA != PSI%NORB'

    ALLOCATE(U(NA,NA))
    U = OVEC

    CALL WFN_INIT(TAU,PSI%NEL,PSI%NORB,PSI%MULT)
    ALLOCATE(T(NA))

    DO J=1,NA
      ! get column vector
      T(:) = U(:,J)
      ! invert T, unitify U(:,J)
      FACT = 1.0D0/T(J)
      DO I=1,NA
        T(I) = -FACT*T(I)
        U(I,J) = 0.0D0
      END DO
      T(J) = FACT
      U(J,J) = 1.0D0
      ! apply T^-1 to U for the next step
      DO K=J+1,NA
        FACT=U(J,K)
        DO I=1,NA
          U(I,K) = U(I,K) + FACT*T(I)
        END DO
        U(J,K) = FACT*T(J)
      END DO
      ! compute |PSI> := (1+E_ij+E_ijmj)|PSI>
      TAU%COEF=(1.5D0-0.5D0*T(J))*PSI%COEF
      DO I=1,NA
        FACT=0.5D0*T(I)
        IF (I.EQ.J) FACT=FACT-0.5D0
        CALL WFN_EX1(FACT,I,J,PSI,TAU)
      END DO
      DO I=1,NA
        FACT=T(I)
        IF (I.EQ.J) FACT=FACT-1.0D0
        CALL WFN_EX1(FACT,I,J,TAU,PSI)
      END DO
    END DO

    DEALLOCATE(T)
    DEALLOCATE(U)
    CALL WFN_FREE(TAU)
  END SUBROUTINE WFN_ORBROT

  SUBROUTINE WFN_ENERGY(NI,NA,NS,D1,D2,ONEINT,TWOINT)
    REAL(REAL64) :: D1(:,:), D2(:,:,:,:)
    REAL(REAL64) :: ONEINT(:,:), TWOINT(:,:,:,:)
    INTEGER :: NI, NA, NS
    INTEGER :: T, U, V, X, I, J
    REAL(REAL64) :: E1, E2, ETOT
    WRITE(*,*) 'Total CASSCF energy:'
    ETOT=-16.704241392292605D0
    WRITE(*,*) '1-el contribution'
    ! inactive
    E1=0.0D0
    DO I=1,NI
      E1=E1+ONEINT(I,I)*2.0D0
    END DO
    WRITE(*,'(1X,A10,F21.14)') 'inactive:', E1
    ETOT=ETOT+E1
    ! active
    E1=0.0D0
    DO T=1,NA
      DO U=1,NA
        IF (D1(T,U).NE.0.0D0) THEN
          E1=E1+ONEINT(NI+T,NI+U)*D1(T,U)
        END IF
      END DO
    END DO
    WRITE(*,'(1X,A10,F21.14)') 'active:', E1
    ETOT=ETOT+E1

    WRITE(*,*) '2-el contribution'
    ! inactive
    E2=0.0D0
    DO I=1,NI
      DO J=1,NI
        E2=E2+2.0D0*TWOINT(I,I,J,J)-TWOINT(I,J,J,I)
      END DO
    END DO
    WRITE(*,'(1X,A10,F21.14)') 'inactive:', E2
    ETOT=ETOT+E2
    ! inactive-active
    E2=0.0D0
    DO I=1,NI
      DO T=1,NA
        DO U=1,NA
          IF (D1(T,U).NE.0.0D0) THEN
            IF (ISNAN(TWOINT(I,I,NI+T,NI+U)).OR.ISNAN(TWOINT(I,NI+U,NI+T,I))) THEN
              WRITE(*,*) 'NaN detected on: ', I, I, NI+T, NI+U
            END IF
            E2=E2+(2.0D0*TWOINT(I,I,NI+T,NI+U)-TWOINT(I,NI+U,NI+T,I))*D1(T,U)
          END IF
        END DO
      END DO
    ENDDO
    WRITE(*,'(1X,A10,F21.14)') 'inact-act:', E2
    ETOT=ETOT+E2
    E2=0.0D0
    DO T=1,NA
      DO U=1,NA
        DO V=1,NA
          DO X=1,NA
            IF (D2(V,T,X,U).NE.0.0D0) THEN
              IF (ISNAN(TWOINT(NI+T,NI+U,NI+V,NI+X))) THEN
                WRITE(*,*) 'NaN detected on: ', NI+T, NI+U, NI+V, NI+X
              END IF
              E2=E2+0.5*TWOINT(NI+T,NI+U,NI+V,NI+X)*D2(V,T,X,U)
            END IF
          END DO
        END DO
      END DO
    ENDDO
    WRITE(*,'(1X,A10,F21.14)') 'active:', E2
    ETOT=ETOT+E2
    WRITE(*,*)
    WRITE(*,'(1X,A10,F21.14)') 'total:', ETOT
    WRITE(*,*)
  END SUBROUTINE WFN_ENERGY

END MODULE WAVEFUNCTION
