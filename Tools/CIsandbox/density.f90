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
MODULE DENSITY
  ! implements 1-electron density matrix
  ! implements 2-electron density matrix

  ! Steven Vancoillie, November 2013, Lund

  USE ISO_FORTRAN_ENV
  USE SECOND_QUANTIZATION
  USE WAVEFUNCTION
  IMPLICIT NONE

CONTAINS

  SUBROUTINE D1_EXC(PSI,D1)
    ! use excitation operator and DDOT to build
    ! the one-electron density matrix
    IMPLICIT NONE
    TYPE(WFN), INTENT(IN) :: PSI
    REAL(REAL64), INTENT(OUT) :: D1(:,:)
    REAL(REAL64), ALLOCATABLE :: SIGMA(:,:)
    INTEGER :: P, Q, IDET, DET, TMP, IRANK
    REAL(REAL64), EXTERNAL :: DDOT

    ALLOCATE(SIGMA(PSI%NDETA,PSI%NDETB))

    DO Q = 1, PSI%NORB
      DO P = 1, PSI%NORB
        SIGMA = 0.0D0

        DET = LEX_INIT(PSI%NELB,PSI%NORB)
        DO IDET = 1, PSI%NDETB
          TMP = EX1(P,Q,DET)
          IF (TMP.NE.-1) THEN
            IRANK=RANK(TMP)
            SIGMA(:,IRANK) = SIGMA(:,IRANK) + PSI%COEF(:,IDET) * SIGN(1,TMP)
          END IF
          DET = LEX_NEXT(DET)
        END DO

        DET = LEX_INIT(PSI%NELA,PSI%NORB)
        DO IDET = 1, PSI%NDETA
          TMP = EX1(P,Q,DET)
          IF (TMP.NE.-1) THEN
            IRANK=RANK(TMP)
            SIGMA(IRANK,:) = SIGMA(IRANK,:) + PSI%COEF(IDET,:) * SIGN(1,TMP)
          END IF
          DET = LEX_NEXT(DET)
        END DO

        D1(P,Q) = DDOT(PSI%NDETA*PSI%NDETB,SIGMA,1,PSI%COEF,1)
      END DO
    END DO
  END SUBROUTINE D1_EXC

  SUBROUTINE D1_ANN(PSI,D1)
    ! use annihilation operators on bra and ket
    ! combined with matrix multiply to build the
    ! one-electron density matrix
    IMPLICIT NONE
    TYPE(WFN), INTENT(IN) :: PSI
    REAL(REAL64), INTENT(OUT) :: D1(:,:)
    REAL(REAL64), ALLOCATABLE :: BRA(:,:,:), KET(:,:,:)
    INTEGER :: P, Q, NORB
    INTEGER :: IA, IB, MDETA, MDETB, NDETSUB
    INTEGER :: DET, TMP, IRANK
    INTEGER :: F
    REAL(REAL64) :: TRACE

    NORB=PSI%NORB

    D1=0.0D0

    ! alpha subspace
    MDETA=BINOM_COEF(PSI%NELA-1,PSI%NORB)
    MDETB=PSI%NDETB
    NDETSUB=MDETA*MDETB

    ALLOCATE(BRA(MDETB,MDETA,NORB))
    ALLOCATE(KET(MDETB,MDETA,NORB))

    KET=0.0D0
    DO P=1,NORB
      DET=LEX_INIT(PSI%NELA,NORB)
      DO IA=1,PSI%NDETA
        TMP=ANN(P,DET)
        IF (TMP.NE.-1) THEN
          IRANK=RANK(TMP)
          F=SIGN(1,TMP)
          DO IB=1,PSI%NDETB
            KET(IB,IRANK,P) = PSI%COEF(IA,IB) * F
          END DO
        END IF
        DET = LEX_NEXT(DET)
      END DO
    END DO
    BRA=KET

    CALL DGEMM('T','N',NORB,NORB,NDETSUB,1.0D0,BRA,NDETSUB,KET,NDETSUB,1.0D0,D1,NORB)

    DEALLOCATE(BRA)
    DEALLOCATE(KET)

    ! beta subspace
    MDETA=PSI%NDETA
    MDETB=BINOM_COEF(PSI%NELB-1,PSI%NORB)
    NDETSUB=MDETA*MDETB

    ALLOCATE(BRA(MDETA,MDETB,NORB))
    ALLOCATE(KET(MDETA,MDETB,NORB))

    KET=0.0D0
    DO P = 1, PSI%NORB
      DET=LEX_INIT(PSI%NELB,PSI%NORB)
      DO IB=1,PSI%NDETB
        TMP=ANN(P,DET)
        IF (TMP.NE.-1) THEN
          IRANK=RANK(TMP)
          F=SIGN(1,TMP)
          DO IA=1,PSI%NDETA
            KET(IA,IRANK,P) = PSI%COEF(IA,IB) * F
          END DO
        END IF
        DET = LEX_NEXT(DET)
      END DO
    END DO
    BRA=KET

    CALL DGEMM('T','N',NORB,NORB,NDETSUB,1.0D0,BRA,NDETSUB,KET,NDETSUB,1.0D0,D1,NORB)

    DEALLOCATE(BRA)
    DEALLOCATE(KET)

  END SUBROUTINE D1_ANN

  SUBROUTINE D2_ANN(PSI,D2)
    ! use annihilation operators on bra and ket
    ! combined with matrix multiply to build the
    ! two-electron density matrix
    !
    ! first compute â_s â_q | psi > for all s,q
    ! then use that as both bra and ket to construct
    ! the matrix D(r,p,s,q) = <psi|E_pqrs|psi>
    ! this is done for each spin combination of
    ! the s and q orbitals and then summed
    IMPLICIT NONE
    TYPE(WFN), INTENT(IN) :: PSI
    REAL(REAL64), INTENT(OUT) :: D2(:,:,:,:)
    REAL(REAL64), ALLOCATABLE :: BRA(:,:,:,:), KET(:,:,:,:)
    INTEGER :: P, Q, NORB, NORB2
    INTEGER :: IA, IB, MDETA, MDETB, NDETSUB
    INTEGER :: DETA, DETB, TMPA, TMPB, IRANKA, IRANKB
    INTEGER :: F, FA, FB
    REAL(REAL64) :: TRACE

    D2 = 0.0D0

    NORB=PSI%NORB
    NORB2=NORB**2

    ! alpha-alpha subspace
    MDETA=BINOM_COEF(PSI%NELA-2,PSI%NORB)
    MDETB=PSI%NDETB
    NDETSUB=MDETA*MDETB

    ALLOCATE(BRA(MDETB,MDETA,NORB,NORB))
    ALLOCATE(KET(MDETB,MDETA,NORB,NORB))

    KET=0.0D0
    DO P=1,NORB
      DO Q=1,NORB
        DETA=LEX_INIT(PSI%NELA,NORB)
        DO IA=1,PSI%NDETA
          TMPA=ANN2(P,Q,DETA)
          IF (TMPA.NE.-1) THEN
            IRANKA=RANK(TMPA)
            FA=SIGN(1,TMPA)
            DO IB=1,PSI%NDETB
              KET(IB,IRANKA,P,Q) = PSI%COEF(IA,IB) * FA
            END DO
          END IF
          DETA = LEX_NEXT(DETA)
        END DO
      END DO
    END DO
    BRA=KET

    CALL DGEMM('T','N',NORB2,NORB2,NDETSUB,1.0D0, &
        & BRA,NDETSUB,KET,NDETSUB,1.0D0,D2,NORB2)

    DEALLOCATE(BRA)
    DEALLOCATE(KET)

    ! beta-beta subspace
    MDETA=PSI%NDETA
    MDETB=BINOM_COEF(PSI%NELB-2,PSI%NORB)
    NDETSUB=MDETA*MDETB
    ALLOCATE(BRA(MDETA,MDETB,NORB,NORB))
    ALLOCATE(KET(MDETA,MDETB,NORB,NORB))

    KET=0.0D0
    DO P=1,NORB
      DO Q=1,NORB
        DETB=LEX_INIT(PSI%NELB,PSI%NORB)
        DO IB=1,PSI%NDETB
          TMPB=ANN2(P,Q,DETB)
          IF (TMPB.NE.-1) THEN
            IRANKB=RANK(TMPB)
            FB=SIGN(1,TMPB)
            DO IA=1,PSI%NDETA
              KET(IA,IRANKB,P,Q) = PSI%COEF(IA,IB) * FB
            END DO
          END IF
          DETB = LEX_NEXT(DETB)
        END DO
      END DO
    END DO
    BRA=KET

    CALL DGEMM('T','N',NORB2,NORB2,NDETSUB,1.0D0, &
        & BRA,NDETSUB,KET,NDETSUB,1.0D0,D2,NORB2)

    DEALLOCATE(BRA)
    DEALLOCATE(KET)

    ! alpha beta
    MDETA=BINOM_COEF(PSI%NELA-1,PSI%NORB)
    MDETB=BINOM_COEF(PSI%NELB-1,PSI%NORB)
    NDETSUB=MDETA*MDETB
    ALLOCATE(BRA(MDETB,MDETA,NORB,NORB))
    ALLOCATE(KET(MDETB,MDETA,NORB,NORB))

    KET=0.0D0
    DO P=1,NORB
      DO Q=1,NORB
        DETA=LEX_INIT(PSI%NELA,NORB)
        DO IA=1,PSI%NDETA
          TMPA=ANN(P,DETA)
          IF (TMPA.NE.-1) THEN
            IRANKA=RANK(TMPA)
            FA=SIGN(1,TMPA)
            DETB=LEX_INIT(PSI%NELB,NORB)
            DO IB=1,PSI%NDETB
              TMPB=ANN(Q,DETB)
              IF (TMPB.NE.-1) THEN
                IRANKB=RANK(TMPB)
                FB=SIGN(1,TMPB)*(-1)**(PSI%NELA)
                KET(IRANKB,IRANKA,P,Q) = PSI%COEF(IA,IB) * FA * FB
              END IF
              DETB = LEX_NEXT(DETB)
            END DO
          END IF
          DETA = LEX_NEXT(DETA)
        END DO
      END DO
    END DO
    BRA=KET

    CALL DGEMM('T','N',NORB2,NORB2,NDETSUB,1.0D0, &
        & BRA,NDETSUB,KET,NDETSUB,1.0D0,D2,NORB2)

    DEALLOCATE(BRA)
    DEALLOCATE(KET)

    ! beta alpha
    MDETA=BINOM_COEF(PSI%NELA-1,PSI%NORB)
    MDETB=BINOM_COEF(PSI%NELB-1,PSI%NORB)
    NDETSUB=MDETA*MDETB
    ALLOCATE(BRA(MDETA,MDETB,NORB,NORB))
    ALLOCATE(KET(MDETA,MDETB,NORB,NORB))

    KET=0.0D0
    DO P=1,NORB
      DO Q=1,NORB
        DETB=LEX_INIT(PSI%NELB,NORB)
        DO IB=1,PSI%NDETB
          TMPB=ANN(P,DETB)
          IF (TMPB.NE.-1) THEN
            IRANKB=RANK(TMPB)
            FB=SIGN(1,TMPB)*(-1)**(PSI%NELA-1)
            DETA=LEX_INIT(PSI%NELA,NORB)
            DO IA=1,PSI%NDETA
              TMPA=ANN(Q,DETA)
              IF (TMPA.NE.-1) THEN
                IRANKA=RANK(TMPA)
                FA=SIGN(1,TMPA)
                KET(IRANKA,IRANKB,P,Q) = PSI%COEF(IA,IB) * FA * FB
              END IF
              DETA = LEX_NEXT(DETA)
            END DO
          END IF
          DETB = LEX_NEXT(DETB)
        END DO
      END DO
    END DO
    BRA=KET

    CALL DGEMM('T','N',NORB2,NORB2,NDETSUB,1.0D0, &
        & BRA,NDETSUB,KET,NDETSUB,1.0D0,D2,NORB2)

    DEALLOCATE(BRA)
    DEALLOCATE(KET)

  END SUBROUTINE D2_ANN

END MODULE DENSITY
