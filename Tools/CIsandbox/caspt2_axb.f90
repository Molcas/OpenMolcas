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
SUBROUTINE CASPT2_AXB(PSI,F,NMO,ONEINT,TWOINT)
  ! Computes the CASPT2 second-order energy of the wavefunction PSI, provided
  ! that the full Fock matrix and necessary 1- and 2-el integrals are given.
  USE ISO_FORTRAN_ENV
  USE SECOND_QUANTIZATION
  USE WAVEFUNCTION
  USE DENSITY
  USE FOCKMATRIX
  USE ORBINT
  IMPLICIT NONE

  TYPE(WFN), INTENT(IN) :: PSI
  TYPE(FMAT), INTENT(IN) :: F
  INTEGER, INTENT(IN) :: NMO
  REAL(REAL64), INTENT(IN) :: ONEINT(NMO,NMO), TWOINT(NMO,NMO,NMO,NMO)

  INTEGER :: NI, NA, NS
  INTEGER :: NACTEL, NTOTEL, MULT

  REAL(REAL64), ALLOCATABLE :: FTOT(:,:)

  TYPE(WFN) :: PSI1, SGM, TAU
  REAL(REAL64), ALLOCATABLE :: D1(:,:), D2(:,:,:,:)
  REAL(REAL64) :: FACT

  REAL(REAL64), ALLOCATABLE :: VPQRS(:,:,:,:), V0(:), X(:)
  REAL(REAL64), ALLOCATABLE :: H0(:,:), S0(:,:), U0(:,:), T0(:,:)
  REAL(REAL64) :: E0, E2
  INTEGER :: ISD, JSD, NSD, MSD

  REAL(REAL64), ALLOCATABLE :: PSI1BRA(:,:,:), PSI1KET(:,:,:), FOCKKET(:,:,:)

  REAL(REAL64), ALLOCATABLE :: WORK(:), EVAL(:)
  INTEGER :: NWORK, INFO

  INTEGER, ALLOCATABLE :: IPIV(:)

  INTEGER :: I, J, P, Q, R, S, IA, IB
  INTEGER :: DETA, DETB, TMPA, TMPB

  REAL(REAL64), EXTERNAL :: DNRM2, DDOT

  ! get appropriate orbital space dimensions from the supplied Fock matrix
  NI=SIZE(F%II,1)
  NA=SIZE(F%AA,1)
  NS=SIZE(F%SS,1)

#ifdef _DEBUG_
  WRITE(*,'(1X,A,I6)') '#inactive  orbitals = ', NI
  WRITE(*,'(1X,A,I6)') '#active    orbitals = ', NA
  WRITE(*,'(1X,A,I6)') '#secondary orbitals = ', NS
  WRITE(*,'(1X,A,I6)') '#molecular orbitals = ', NMO
#endif

  ! some sanity checks
  IF (PSI%NORB.NE.NA) STOP 'CASPT2_AXB: PSI%NORB does not match #active orbitals'
  IF (NMO.NE.NI+NA+NS) STOP 'CASPT2_AXB: supplied total #orbitals does not compute'

  NACTEL=PSI%NEL
  NTOTEL=2*NI+NACTEL
  MULT=PSI%MULT

  ! put Fock matrix in general matrix format
  ALLOCATE(FTOT(NMO,NMO))
  CALL FMAT_COPY(F,FTOT,NMO)

  ! first, set up PSI1 as a (sparse) full-CI wavefunction for which we will use
  ! projectors to handle the restriction of single-double replacement states
  CALL WFN_INIT(PSI1,NTOTEL,NMO,MULT)

  PSI1%COEF=0.0D0
  CALL WFN_COPY0(NI,PSI,PSI1)

#ifdef _DEBUG_
  WRITE(*,'(1X,A)') 'The 0th-order wavefunction in FCI space:'
  CALL WFN_PRINT(PSI1,1.0D-15)
#endif

  ! compute the CASSCF energy for the 0-th order wavefunction
  ALLOCATE(D1(NMO,NMO))
  ALLOCATE(D2(NMO,NMO,NMO,NMO))
  CALL D1_ANN(PSI1,D1)
  CALL D2_ANN(PSI1,D2)
  
  CALL WFN_ENERGY(0,NMO,0,D1,D2,ONEINT,TWOINT)

  CALL WFN_INIT(SGM,NTOTEL,NMO,MULT)
  CALL WFN_INIT(TAU,NTOTEL,NMO,MULT)

  ! set up the equation system Ax=b as (H0 -E0*S0) X = -V0 with:
  ! H0_ij = <i|Ĥ_0|j>
  ! E0    = <0|Ĥ_0|0>
  ! S0_ij = <i|j>
  ! V0_i  = <i|Ĥ|0>
  ! |i>, |j> = P_SD Ê_pqrs |0>

  ! First, compute the 0-th order energy and construct the right hand side
  ! vector -V0, which is computed for each |i> that is non-zero. This still
  ! leaves an overcomplete basis for the VSD space, but that is solved later by
  ! diagonalization of S0 and removal of any linear dependency before solving
  ! Ax=b.

  ! step 1: |TAU> = Ĥ|0>; |SGM> = Ĥ_0|0>
  TAU%COEF=0.0D0
  SGM%COEF=0.0D0

  DETB=LEX_INIT(PSI%NELB,PSI%NORB)
  DO IB=1,PSI%NDETB
    TMPB=ISHFT(DETB,NI)+(2**NI-1)
    DETA=LEX_INIT(PSI%NELA,PSI%NORB)
    DO IA=1,PSI%NDETA
      TMPA=ISHFT(DETA,NI)+(2**NI-1)

      FACT=PSI%COEF(IA,IB)
      DO P=1,NMO
        DO Q=1,NMO
          ! h_pq, g_pqrs for TAU
          CALL DET_EX1(FACT*ONEINT(P,Q),P,Q,TMPA,TMPB,TAU)
          DO R=1,NMO
            DO S=1,NMO
              CALL DET_EX2(FACT*0.5D0*TWOINT(P,Q,R,S),P,Q,R,S,TMPA,TMPB,TAU)
            END DO
          END DO
          ! f_pq for SGM
          CALL DET_EX1(FACT*FTOT(P,Q),P,Q,TMPA,TMPB,SGM)
        END DO
      END DO

      DETA=LEX_NEXT(DETA)
    END DO
    DETB=LEX_NEXT(DETB)
  END DO

#ifdef _DEBUG_
  WRITE(*,'(1X,A)') 'wavefunction Ĥ|0>'
  CALL WFN_PRINT(TAU,1.0D-6)
  WRITE(*,'(1X,A)') 'wavefunction Ĥ_0|0>'
  CALL WFN_PRINT(SGM,1.0D-6)
#endif

  ! step 2: compute 0-th order energy E0 = <0|SGM>
  E0 = 0.0D0

  DETB=LEX_INIT(PSI%NELB,PSI%NORB)
  DO IB=1,PSI%NDETB
    TMPB=ISHFT(DETB,NI)+(2**NI-1)
    DETA=LEX_INIT(PSI%NELA,PSI%NORB)
    DO IA=1,PSI%NDETA
      TMPA=ISHFT(DETA,NI)+(2**NI-1)
  
      E0 = E0 + PSI%COEF(IA,IB) * SGM%COEF(RANK_(TMPA),RANK_(TMPB))

      DETA=LEX_NEXT(DETA)
    END DO
    DETB=LEX_NEXT(DETB)
  END DO

  WRITE(*,'(1X,A32,F21.14)') "0-th order energy E0 = ", E0

  ! step 3: construct V0 = <i|TAU>, store the computed functions |i> in PSI1KET
  ! and PSI1BRA
  ALLOCATE(VPQRS(NMO,NMO,NMO,NMO))

  NSD=0
  VPQRS=0.0D0
  DO S=1,NMO
    DO R=1,NMO
      DO Q=1,NMO
        DO P=1,NMO

          ! determine |i> = E_pqrs |Psi0>
          SGM%COEF=0.0D0
          DETB=LEX_INIT(PSI%NELB,PSI%NORB)
          DO IB=1,PSI%NDETB
            TMPB=ISHFT(DETB,NI)+(2**NI-1)
            DETA=LEX_INIT(PSI%NELA,PSI%NORB)
            DO IA=1,PSI%NDETA
              TMPA=ISHFT(DETA,NI)+(2**NI-1)

              CALL DET_EX2(PSI%COEF(IA,IB),P,Q,R,S,TMPA,TMPB,SGM)

              DETA=LEX_NEXT(DETA)
            END DO
            DETB=LEX_NEXT(DETB)
          END DO

          PSI1%COEF=0.0D0
          CALL WFN_PSD(NI,NA,NS,SGM,PSI1)

          IF (DNRM2(PSI1%NDET,PSI1%COEF,1).NE.0.0D0) THEN
            VPQRS(P,Q,R,S)=DDOT(PSI1%NDET,PSI1%COEF,1,TAU%COEF,1)
            NSD=NSD+1
          END IF

        END DO
      END DO
    END DO
  END DO

  ALLOCATE(V0(NSD))

  ISD=0
  DO S=1,NMO
    DO R=1,NMO
      DO Q=1,NMO
        DO P=1,NMO
          IF (ABS(VPQRS(P,Q,R,S)).GT.0.0D0) THEN
            ISD=ISD+1
            V0(ISD)=VPQRS(P,Q,R,S)
          END IF
        END DO
      END DO
    END DO
  END DO

  ! Now compute H0 and S0 matrices

  ALLOCATE(H0(NSD,NSD))
  ALLOCATE(S0(NSD,NSD))

  H0 = 0.0D0
  S0 = 0.0D0

  ALLOCATE(PSI1BRA(NSD,PSI1%NDETA,PSI1%NDETB))
  ALLOCATE(PSI1KET(PSI1%NDETA,PSI1%NDETB,NSD))
  ALLOCATE(FOCKKET(PSI1%NDETA,PSI1%NDETB,NSD))

  ISD=0
  DO S=1,NMO
    DO R=1,NMO
      DO Q=1,NMO
        DO P=1,NMO

          SGM%COEF=0.0D0
          ! determine |i> = E_pqrs |Psi0>
          DETB=LEX_INIT(PSI%NELB,PSI%NORB)
          DO IB=1,PSI%NDETB
            TMPB=ISHFT(DETB,NI)+(2**NI-1)
            DETA=LEX_INIT(PSI%NELA,PSI%NORB)
            DO IA=1,PSI%NDETA
              TMPA=ISHFT(DETA,NI)+(2**NI-1)

              CALL DET_EX2(PSI%COEF(IA,IB),P,Q,R,S,TMPA,TMPB,SGM)

              DETA=LEX_NEXT(DETA)
            END DO
            DETB=LEX_NEXT(DETB)
          END DO
          PSI1%COEF=0.0D0
          CALL WFN_PSD(NI,NA,NS,SGM,PSI1)

          IF (DNRM2(PSI1%NDET,PSI1%COEF,1).NE.0.0D0) THEN
            ISD=ISD+1
            ! put |i> = E_pqrs |0> into PSI1BRA and PSI1KET
            PSI1BRA(ISD,:,:) = PSI1%COEF
            PSI1KET(:,:,ISD) = PSI1%COEF

            ! compute Ĥ_0|i> = Sum_ij f_ij Ê_ij |i>
            SGM%COEF=0.0D0
            DETB=LEX_INIT(PSI1%NELB,PSI1%NORB)
            DO IB=1,PSI1%NDETB
              DETA=LEX_INIT(PSI1%NELA,PSI1%NORB)
              DO IA=1,PSI1%NDETA

                DO J=1,NMO
                  DO I=1,NMO
                    ! accumulate F_ij Ê_ij |i>
                    CALL DET_EX1(PSI1%COEF(IA,IB)*FTOT(I,J),I,J,DETA,DETB,SGM)
                  END DO
                END DO

                DETA=LEX_NEXT(DETA)
              END DO
              DETB=LEX_NEXT(DETB)
            END DO
            ! put Ĥ_0|i> into FOCKKET
            FOCKKET(:,:,ISD) = SGM%COEF
          END IF

        END DO
      END DO
    END DO
  END DO

  CALL DGEMM('N','N',NSD,NSD,PSI1%NDET,1.0D0, &
       & PSI1BRA,NSD,PSI1KET,PSI1%NDET, &
       & 0.0D0,S0,NSD)

  CALL DGEMM('N','N',NSD,NSD,PSI1%NDET,1.0D0, &
       & PSI1BRA,NSD,FOCKKET,PSI1%NDET, &
       & 0.0D0,H0,NSD)

  DEALLOCATE(PSI1BRA)
  DEALLOCATE(PSI1KET)
  DEALLOCATE(FOCKKET)

  ! remove linear dependencies by diagonalizing the overlap matrix, construct a
  ! transformation matrix to orthogonal basis, and transform S0, H0 and V0
  ! accordingly

  ALLOCATE(U0(NSD,NSD))
  U0=S0

  NWORK=NSD**2
  ALLOCATE(WORK(NWORK))
  ALLOCATE(EVAL(NSD))

  call dsyev_('V','U',NSD,U0,NSD,EVAL,WORK,NWORK,INFO)
  IF (INFO.NE.0) STOP 'CASPT2_UDL: diagonalization of S0 failed'

  ! go through columns of U0, keep non-zero eigenvalues
  MSD=0
  DO ISD=1,NSD
    IF (EVAL(ISD).GT.1.0D-6) THEN
      MSD=MSD+1
    END IF
  END DO
  WRITE(*,'(1X,A,I6)') '#non-zero overlap eigenvalues = ', MSD

  ALLOCATE(T0(NSD,MSD))
  ISD=0
  DO JSD=1,NSD
    IF (EVAL(JSD).GT.1.0D-6) THEN
      ISD=ISD+1
      T0(:,ISD) = U0(:,JSD)/SQRT(EVAL(JSD))
    END IF
  END DO

  ! transform equation system
  H0(1:MSD,1:MSD) = MATMUL(TRANSPOSE(T0),MATMUL(H0,T0))
  S0(1:MSD,1:MSD) = MATMUL(TRANSPOSE(T0),MATMUL(S0,T0))

  ! Solve (H0 - E0 S0) X = -V0

  H0 = H0 - E0 * S0

  ALLOCATE(IPIV(NSD))

  ! solution vector is initially filled with the RHS = -V0
  ALLOCATE(X(NSD))
  X(1:MSD) = -1.0D0 * MATMUL(TRANSPOSE(T0),V0)

  call dgesv_(MSD,1,H0,NSD,IPIV,X,NSD,INFO)
  IF (INFO.NE.0) THEN
    WRITE(*,*) 'INFO = ', INFO
    STOP 'Failed to solve linear equation system'
  END IF

  X=MATMUL(T0,X(1:MSD))

  ! compute the energy

  ! method 1: E2 = <V0|X>
  E2 = DDOT(NSD,V0,1,X,1)
  WRITE(*,'(1X,A32,F21.14)') 'E2 as <V0|X> = ', E2

  ! method 2: E2 = <0|Ĥ|Psi1>
  ! first, form the first-order wavefunction |Psi1>
  ISD=0
  SGM%COEF=0.0D0
  DO S=1,NMO
    DO R=1,NMO
      DO Q=1,NMO
        DO P=1,NMO

          IF (VPQRS(P,Q,R,S).NE.0.0D0) THEN
            ISD=ISD+1

            DETB=LEX_INIT(PSI%NELB,PSI%NORB)
            DO IB=1,PSI%NDETB
              TMPB=ISHFT(DETB,NI)+(2**NI-1)
              DETA=LEX_INIT(PSI%NELA,PSI%NORB)
              DO IA=1,PSI%NDETA
                TMPA=ISHFT(DETA,NI)+(2**NI-1)

                CALL DET_EX2(PSI%COEF(IA,IB)*X(ISD),P,Q,R,S,TMPA,TMPB,SGM)

                DETA=LEX_NEXT(DETA)
              END DO
              DETB=LEX_NEXT(DETB)
            END DO
            
          END IF

        END DO
      END DO
    END DO
  END DO
  PSI1%COEF=0.0D0
  CALL WFN_PSD(NI,NA,NS,SGM,PSI1)

  !WRITE(*,*) 'Final first-order wavefunction from Ax=b method:'
  !CALL WFN_PRINT(PSI1,1.0D-8)

  ! TAU still contains Ĥ|0>, so contract with |Psi1>
  E2 = DDOT(PSI1%NDET,PSI1%COEF,1,TAU%COEF,1)
  WRITE(*,'(1X,A32,F21.14)') 'E2 as <Psi1|H|0> = ', E2

END SUBROUTINE
