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
SUBROUTINE CASPT2_UDL(PSI,F,NAO,NMO,CMO,ONEINT,TWOINT)
  ! Computes the CASPT2 second-order energy of the wavefunction PSI, provided
  ! that the full Fock matrix and necessary 1- and 2-el integrals are given.

  ! Uses UDL decomposition of the Fock matrix and transforms the RHS according
  ! to the upper and lower triangular blocks Q of U and L by means of the
  ! exponential operator exp(Q).

  USE ISO_FORTRAN_ENV
  USE SECOND_QUANTIZATION
  USE WAVEFUNCTION
  USE DENSITY
  USE FOCKMATRIX
  USE ORBINT
  IMPLICIT NONE

  TYPE(WFN), INTENT(INOUT) :: PSI
  TYPE(FMAT), INTENT(INOUT) :: F
  INTEGER, INTENT(IN) :: NAO, NMO
  REAL(REAL64), INTENT(INOUT) :: ONEINT(NMO,NMO), TWOINT(NMO,NMO,NMO,NMO)
  REAL(REAL64), INTENT(INOUT) :: CMO(NAO,NMO)

  INTEGER :: NI, NA, NS
  INTEGER :: NACTEL, NTOTEL, MULT

  TYPE(WFN) :: PSI1, SGM, TAU, RHS
  REAL(REAL64), ALLOCATABLE :: D1(:,:), D2(:,:,:,:), FTOT(:,:)
  REAL(REAL64) :: FACT

  REAL(REAL64), ALLOCATABLE :: VPQRS(:,:,:,:), V0(:), X(:)
  REAL(REAL64), ALLOCATABLE :: H0(:,:), S0(:,:), U0(:,:), T0(:,:)
  REAL(REAL64) :: E0, E2
  INTEGER :: ISD, JSD, NSD, MSD

  REAL(REAL64), ALLOCATABLE :: PSI1BRA(:,:,:), PSI1KET(:,:,:), FOCKKET(:,:,:)

  REAL(REAL64), ALLOCATABLE :: WORK(:), EVAL(:)
  INTEGER :: NWORK, INFO

  INTEGER, ALLOCATABLE :: IPIV(:)

  INTEGER :: A, B, I, J, T, U, P, Q, R, S, IA, IB
  INTEGER :: DETA, DETB, TMPA, TMPB
  INTEGER, ALLOCATABLE :: IDXA(:), IDXB(:), ADET(:), BDET(:)

  REAL(REAL64), EXTERNAL :: DNRM2, DDOT

  INTEGER :: ICNT

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

  ! reconstruct density matrices
  ALLOCATE(D1(NA,NA))
  ALLOCATE(D2(NA,NA,NA,NA))
  CALL D1_ANN(PSI,D1)
  CALL D2_ANN(PSI,D2)

  WRITE(*,*) 'Before inactive/active orbital rotation:'
  WRITE(*,*)
  CALL WFN_PRINT(PSI,0.0D0)
  CALL WFN_ENERGY(NI,NA,NS,D1,D2,ONEINT,TWOINT)

  ! copy Fock matrix and apply UDL decomposition after shifting by E0!
  E0 = -4.39765297778118
  DO I=1,NI
    F%II(I,I)=F%II(I,I)-E0/NTOTEL
  END DO
  DO T=1,NA
    F%AA(T,T)=F%AA(T,T)-E0/NTOTEL
  END DO
  DO A=1,NS
    F%SS(A,A)=F%SS(A,A)-E0/NTOTEL
  END DO

  CALL FMAT_UDL(F)
  WRITE(*,*) 'After UDL decomposition'
  CALL FMAT_PRINT(F)

  ! store diagonal in FTOT
  ALLOCATE(FTOT(NMO,NMO))
  FTOT=0.0D0
  DO I=1,NI
    FTOT(I,I)=F%DI(I)
  END DO
  DO T=1,NA
    FTOT(NI+T,NI+T)=F%DA(T)
  END DO
  DO A=1,NS
    FTOT(NI+NA+A,NI+NA+A)=F%DS(A)
  END DO

!#ifdef _DEBUG_
  DO P=1,NMO
    WRITE(*,'(1X,7(2X,F21.14))') FTOT(P,:)
  END DO
!#endif

  ! copy and transform the wavefunction accordingly
  CALL WFN_ORBROT(PSI,F%AA)
  CALL ORBINT_TRANSFORM(CMO,ONEINT,TWOINT,F%II,F%AA,F%SS)

  ! reconstruct density matrices
  CALL D1_ANN(PSI,D1)
  CALL D2_ANN(PSI,D2)

  WRITE(*,*) 'After inactive/active orbital rotation:'
  WRITE(*,*)
  CALL WFN_PRINT(PSI,0.0D0)
  CALL WFN_ENERGY(NI,NA,NS,D1,D2,ONEINT,TWOINT)

  ! first, set up PSI1 as a (sparse) full-CI wavefunction for which we will use
  ! projectors to handle the restriction of single-double replacement states
  CALL WFN_INIT(PSI1,NTOTEL,NMO,MULT)

  PSI1%COEF=0.0D0
  CALL WFN_COPY0(NI,PSI,PSI1)

#ifdef _DEBUG_
  WRITE(*,'(1X,A)') 'The 0th-order wavefunction in FCI space:'
  CALL WFN_PRINT(PSI1,1.0D-15)
#endif

  DEALLOCATE(D1)
  DEALLOCATE(D2)

  ! compute the CASSCF energy for the 0-th order wavefunction
  ALLOCATE(D1(NMO,NMO))
  ALLOCATE(D2(NMO,NMO,NMO,NMO))
  CALL D1_ANN(PSI1,D1)
  CALL D2_ANN(PSI1,D2)
  
  CALL WFN_ENERGY(0,NMO,0,D1,D2,ONEINT,TWOINT)

  CALL WFN_INIT(SGM,NTOTEL,NMO,MULT)
  CALL WFN_INIT(TAU,NTOTEL,NMO,MULT)
  CALL WFN_INIT(RHS,NTOTEL,NMO,MULT)

  ! set up the RHS -V0 with V0_i  = <i|Ĥ|0>

  ! First, compute the 0-th order energy and construct the right hand side
  ! vector -V0, which is computed for each |i> that is non-zero. This still
  ! leaves an overcomplete basis for the VSD space, but that is solved later by
  ! diagonalization of S0 and removal of any linear dependency before solving
  ! Ax=b.

  ! step 1: |TAU> = Ĥ|0>
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
        END DO
      END DO

      DETA=LEX_NEXT(DETA)
    END DO
    DETB=LEX_NEXT(DETB)
  END DO

#ifdef _DEBUG_
  WRITE(*,'(1X,A)') 'wavefunction Ĥ|0>'
  CALL WFN_PRINT(TAU,1.0D-6)
#endif

  ! step 2: construct V0 = <i|TAU>, store the computed functions |i> in PSI1KET
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

          IF (DNRM2_(PSI1%NDET,PSI1%COEF,1).NE.0.0D0) THEN
            VPQRS(P,Q,R,S)=DDOT_(PSI1%NDET,PSI1%COEF,1,TAU%COEF,1)
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

  DEALLOCATE(VPQRS)

  ALLOCATE(H0(NSD,NSD))
  ALLOCATE(S0(NSD,NSD))

  H0 = 0.0D0
  S0 = 0.0D0

  ! construct E_pqrs|0> in PSI1KET, PSI1BRA and

  ALLOCATE(PSI1BRA(NSD,PSI1%NDETA,PSI1%NDETB))
  ALLOCATE(PSI1KET(PSI1%NDETA,PSI1%NDETB,NSD))
  ALLOCATE(FOCKKET(PSI1%NDETA,PSI1%NDETB,NSD))

  ISD=0
  DO S=1,NMO
    DO R=1,NMO
      DO Q=1,NMO
        DO P=1,NMO

          SGM%COEF=0.0D0
          ! determine |i> = E_pqrs |0>
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

          IF (DNRM2_(PSI1%NDET,PSI1%COEF,1).NE.0.0D0) THEN
            ISD=ISD+1
            ! put |i> = E_pqrs |0> into PSI1BRA and PSI1KET
            PSI1BRA(ISD,:,:) = PSI1%COEF
            PSI1KET(:,:,ISD) = PSI1%COEF

            ! compute Ĥ_0|i> = Sum_pp f_pp Ê_pp |i>
            SGM%COEF=0.0D0
            DETB=LEX_INIT(PSI1%NELB,PSI1%NORB)
            DO IB=1,PSI1%NDETB
              DETA=LEX_INIT(PSI1%NELA,PSI1%NORB)
              DO IA=1,PSI1%NDETA

                DO I=1,NMO
                  ! accumulate F_ii Ê_ii |i>
                  CALL DET_EX1(PSI1%COEF(IA,IB)*FTOT(I,I),I,I,DETA,DETB,SGM)
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

  ! preparation: determine T^T S T = 1, to remove linear dependencies by
  ! diagonalizing the overlap matrix, construct a transformation matrix to
  ! orthogonal basis, and transform X, S0, H0 and V0 accordingly

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

  DEALLOCATE(U0)

  ALLOCATE(X(NSD))
  X=-V0

  ! compute RHS = Sum_i X_i |i>
  RHS%COEF=0.0D0
  DO ISD=1,NSD
    RHS%COEF=RHS%COEF+X(ISD)*PSI1KET(:,:,ISD)
  END DO

  ! apply first series of operators exp(-Q_IA) exp(-Q_AS) exp(-Q_IS)

  ! exp(-Q_IS) = 1 - Sum_IS Q_IS E_IS + 0.5 ...
  SGM%COEF=0.0D0
  DO A=1,NS
    DO I=1,NI
      CALL WFN_EX1(-F%IS(I,A),I,NI+NA+A,RHS,SGM)
    END DO
  END DO
  PSI1%COEF=SGM%COEF
  DO A=1,NS
    DO I=1,NI
      CALL WFN_EX1(-0.5D0*F%IS(I,A),I,NI+NA+A,PSI1,SGM)
    END DO
  END DO
  RHS%COEF=RHS%COEF+SGM%COEF

  ! exp(-Q_AS) = 1 - Sum_AS Q_AS E_AS + ...
  SGM%COEF=0.0D0
  DO A=1,NS
    DO T=1,NA
      CALL WFN_EX1(-F%AS(T,A),NI+T,NI+NA+A,RHS,SGM)
    END DO
  END DO
  PSI1%COEF=SGM%COEF
  DO A=1,NS
    DO T=1,NA
      CALL WFN_EX1(-0.5D0*F%AS(T,A),NI+T,NI+NA+A,PSI1,SGM)
    END DO
  END DO
  RHS%COEF=RHS%COEF+SGM%COEF

  ! exp(-Q_IA) = 1 - Sum_IA Q_IA E_IA + ...
  SGM%COEF=0.0D0
  DO T=1,NA
    DO I=1,NI
      CALL WFN_EX1(-F%IA(I,T),I,NI+T,RHS,SGM)
    END DO
  END DO
  PSI1%COEF=0.0D0
  CALL WFN_PSD(NI,NA,NS,SGM,PSI1)
  !PSI1%COEF=SGM%COEF
  DO T=1,NA
    DO I=1,NI
      CALL WFN_EX1(-0.5D0*F%IA(I,T),I,NI+T,PSI1,SGM)
    END DO
  END DO
  RHS%COEF=RHS%COEF+SGM%COEF

  ! contract into X_i = <i|exp(-Q_IA)|RHS> and
  ! transform to adapt for the overcomplete basis
  CALL DGEMV_('N',NSD,PSI1%NDET,1.0D0, &
       & PSI1BRA,NSD,RHS%COEF,1,0.0D0,X,1)
  X = MATMUL(MATMUL(T0,TRANSPOSE(T0)),X)


  ! apply diagonal inverse, e.g. solve H0' XNEW = XOLD

  ! transform equation system
  H0(1:MSD,1:MSD) = MATMUL(TRANSPOSE(T0),MATMUL(H0,T0))

  ! solution vector is already filled with RHS, just transform it
  X(1:MSD) = MATMUL(TRANSPOSE(T0),X)

  ALLOCATE(IPIV(NSD))

  call dgesv_(MSD,1,H0,NSD,IPIV,X,NSD,INFO)
  IF (INFO.NE.0) THEN
    WRITE(*,*) 'INFO = ', INFO
    STOP 'Failed to solve linear equation system'
  END IF

  DEALLOCATE(IPIV)

  X=MATMUL(T0,X(1:MSD))

  ! Now, continue with the transformation steps for L

  ! compute new RHS = Sum_i X_i |i>
  RHS%COEF=0.0D0
  DO ISD=1,NSD
    RHS%COEF=RHS%COEF+X(ISD)*PSI1KET(:,:,ISD)
  END DO

  ! apply second series of operators exp(-Q_SI) exp(-Q_SA) exp(-Q_SI)

  ! exp(-Q_AI) = 1 - Sum_AI Q_AI E_AI
  SGM%COEF=0.0D0
  DO I=1,NI
    DO T=1,NA
      CALL WFN_EX1(-F%AI(T,I),NI+T,I,RHS,SGM)
    END DO
  END DO
  PSI1%COEF=SGM%COEF
  DO I=1,NI
    DO T=1,NA
      CALL WFN_EX1(-0.5D0*F%AI(T,I),NI+T,I,PSI1,SGM)
    END DO
  END DO
  RHS%COEF=RHS%COEF+SGM%COEF

  ! exp(-Q_SA) = 1 - Sum_SA Q_SA E_SA + ...
  SGM%COEF=0.0D0
  DO T=1,NA
    DO A=1,NS
      CALL WFN_EX1(-F%SA(A,T),NI+NA+A,NI+T,RHS,SGM)
    END DO
  END DO
  PSI1%COEF=SGM%COEF
  DO T=1,NA
    DO A=1,NS
      CALL WFN_EX1(-0.5D0*F%SA(A,T),NI+NA+A,NI+T,PSI1,SGM)
    END DO
  END DO
  RHS%COEF=RHS%COEF+SGM%COEF

  ! exp(-Q_SI) = 1 - Sum_SI Q_SI E_SI + ...
  SGM%COEF=0.0D0
  DO I=1,NI
    DO A=1,NS
      CALL WFN_EX1(-F%SI(A,I),NI+NA+A,I,RHS,SGM)
    END DO
  END DO
  PSI1%COEF=SGM%COEF
  DO I=1,NI
    DO A=1,NS
      CALL WFN_EX1(-0.5D0*F%SI(A,I),NI+NA+A,I,PSI1,SGM)
    END DO
  END DO
  RHS%COEF=RHS%COEF+SGM%COEF

  ! contract into X_i = <i|exp(-Q_IA)|RHS> and
  ! transform to adapt for the overcomplete basis
  CALL DGEMV_('N',NSD,PSI1%NDET,1.0D0, &
       & PSI1BRA,NSD,RHS%COEF,1,0.0D0,X,1)
  X = MATMUL(MATMUL(T0,TRANSPOSE(T0)),X)



  ! compute the energy

  ! method 1: E2 = <V0|X>
  E2 = DDOT_(NSD,V0,1,X,1)
  WRITE(*,'(1X,A32,F21.14)') 'E2 as <V0|X> = ', E2

  ! method 2: E2 = <0|Ĥ|Psi1>
  ! first, form the first-order wavefunction |Psi1>
  PSI1%COEF=0.0D0
  DO ISD=1,NSD
    PSI1%COEF=PSI1%COEF+X(ISD)*PSI1KET(:,:,ISD)
  END DO

  !WRITE(*,*) 'Final first-order wavefunction from UDL method:'
  !CALL WFN_PRINT(PSI1,1.0D-8)

  ! TAU still contains Ĥ|0>, so contract with |Psi1>
  E2 = DDOT_(PSI1%NDET,PSI1%COEF,1,TAU%COEF,1)
  WRITE(*,'(1X,A32,F21.14)') 'E2 as <Psi1|H|0> = ', E2

  DEALLOCATE(H0)
  DEALLOCATE(S0)
  DEALLOCATE(T0)
  DEALLOCATE(V0)
  DEALLOCATE(X)
  DEALLOCATE(PSI1BRA)
  DEALLOCATE(PSI1KET)
  DEALLOCATE(FOCKKET)



  ! Same UDL decomposition procedure, but in a determinant basis instead of an
  ! excitation parameter space. This has the advantage of being simpler, faster,
  ! and the fact that there is an orthonormal basis in place already.

  ! step 1: |TAU> = Ĥ|0>
  TAU%COEF=0.0D0
  SGM%COEF=0.0D0
  RHS%COEF=0.0D0

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
        END DO
      END DO

      DETA=LEX_NEXT(DETA)
    END DO
    DETB=LEX_NEXT(DETB)
  END DO

  ! project out the non-SD determinants from H|0>
  CALL WFN_PSD(NI,NA,NS,TAU,RHS)

#ifdef _DEBUG_
  WRITE(*,'(1X,A)') 'wavefunction Ĥ|0>'
  CALL WFN_PRINT(TAU,1.0D-6)
#endif

  ! construct PSI1 = Sum_i |i>, store valid indices in separate index tables
  PSI1%COEF=0.0D0
  SGM%COEF=1.0D0
  CALL WFN_PSD(NI,NA,NS,SGM,PSI1)

  NSD=0
  DO IB=1,PSI1%NDETB
    DO IA=1,PSI1%NDETA
      IF (PSI1%COEF(IA,IB).NE.0.0D0) NSD=NSD+1
    END DO
  END DO

  WRITE(*,'(1X,A32,I6)') '#SD determinants in Psi1 = ', NSD

  ! store indices and determinants of non-zero Psi1 contributions in separate
  ! arrays
  ALLOCATE(IDXA(NSD))
  ALLOCATE(IDXB(NSD))
  ALLOCATE(ADET(NSD))
  ALLOCATE(BDET(NSD))

  ISD=0
  TMPB=LEX_INIT(PSI1%NELB,PSI1%NORB)
  DO IB=1,PSI1%NDETB
    TMPA=LEX_INIT(PSI1%NELA,PSI1%NORB)
    DO IA=1,PSI1%NDETA
      IF (PSI1%COEF(IA,IB).NE.0.0D0) THEN
        ISD=ISD+1
        IDXA(ISD)=IA
        IDXB(ISD)=IB
        ADET(ISD)=TMPA
        BDET(ISD)=TMPB
      END IF
      TMPA=LEX_NEXT(TMPA)
    END DO
    TMPB=LEX_NEXT(TMPB)
  END DO

  ! step 2: construct V0, V0_i = <i|H|0>

  ALLOCATE(V0(NSD))

  DO ISD=1,NSD
    IA=IDXA(ISD)
    IB=IDXB(ISD)
    V0(ISD)=RHS%COEF(IA,IB)
  END DO

  ALLOCATE(H0(NSD,NSD))
  ALLOCATE(S0(NSD,NSD))

  H0 = 0.0D0
  S0 = 0.0D0

  ! construct elements <i|H_0|j>
  DO JSD=1,NSD
    TMPA=ADET(JSD)
    TMPB=BDET(JSD)
    SGM%COEF=0.0D0
    DO I=1,NMO
      ! accumulate Sum_i F_ii Ê_ii |j>
      CALL DET_EX1(FTOT(I,I),I,I,TMPA,TMPB,SGM)
    END DO
    DO ISD=1,NSD
      IA=IDXA(ISD)
      IB=IDXB(ISD)
      H0(ISD,JSD)=SGM%COEF(IA,IB)
    END DO
  END DO

  ! set RHS to -<i|H|0>
  RHS%COEF=-RHS%COEF

  ! series of deexcitation operators exp(-Q_IA) exp(-Q_AS) exp(-Q_IS)

  ! exp(-Q_IS) = 1 - Sum_IS Q_IS E_IS + 0.5 ...
  SGM%COEF=0.0D0
  DO A=1,NS
    DO I=1,NI
      CALL WFN_EX1(-F%IS(I,A),I,NI+NA+A,RHS,SGM)
    END DO
  END DO
  PSI1%COEF=SGM%COEF
  DO A=1,NS
    DO I=1,NI
      CALL WFN_EX1(-0.5D0*F%IS(I,A),I,NI+NA+A,PSI1,SGM)
    END DO
  END DO
  RHS%COEF=RHS%COEF+SGM%COEF

  ! exp(-Q_AS) = 1 - Sum_AS Q_AS E_AS + ...
  SGM%COEF=0.0D0
  DO A=1,NS
    DO T=1,NA
      CALL WFN_EX1(-F%AS(T,A),NI+T,NI+NA+A,RHS,SGM)
    END DO
  END DO
  PSI1%COEF=SGM%COEF
  DO A=1,NS
    DO T=1,NA
      CALL WFN_EX1(-0.5D0*F%AS(T,A),NI+T,NI+NA+A,PSI1,SGM)
    END DO
  END DO
  RHS%COEF=RHS%COEF+SGM%COEF

  ! exp(-Q_IA) = 1 - Sum_IA Q_IA E_IA + ...
  SGM%COEF=0.0D0
  DO T=1,NA
    DO I=1,NI
      CALL WFN_EX1(-F%IA(I,T),I,NI+T,RHS,SGM)
    END DO
  END DO
  PSI1%COEF=SGM%COEF
  DO T=1,NA
    DO I=1,NI
      CALL WFN_EX1(-0.5D0*F%IA(I,T),I,NI+T,PSI1,SGM)
    END DO
  END DO
  RHS%COEF=RHS%COEF+SGM%COEF

  ALLOCATE(X(NSD))

  ! contract into X_i = <i|exp(-Q_IA)|RHS>
  DO ISD=1,NSD
    IA=IDXA(ISD)
    IB=IDXB(ISD)
    X(ISD)=RHS%COEF(IA,IB)
  END DO

  ! diagonal resolvent step

  ! method 1: invert H0 and multiply

  !ALLOCATE(U0(NSD,NSD))
  !U0=H0

  !NWORK=NSD**2
  !ALLOCATE(WORK(NWORK))
  !ALLOCATE(EVAL(NSD))

  !call dsyev_('V','U',NSD,U0,NSD,EVAL,WORK,NWORK,INFO)
  !IF (INFO.NE.0) STOP 'CASPT2_UDL: diagonalization of S0 failed'

  !S0=0.0D0
  !DO ISD=1,NSD
  !  S0(ISD,ISD)=1/EVAL(ISD)
  !END DO
  !H0=MATMUL(MATMUL(U0,S0),TRANSPOSE(U0))
  !X=MATMUL(H0,X)

  ! method 2: solve linear system H0 C = X

  ALLOCATE(IPIV(NSD))
  call dgesv_(NSD,1,H0,NSD,IPIV,X,NSD,INFO)
  IF (INFO.NE.0) THEN
    WRITE(*,*) 'INFO = ', INFO
    STOP 'Failed to solve linear equation system'
  END IF

  ! Now, continue with the transformation steps for L

  ! compute new RHS = Sum_i X_i |i>
  RHS%COEF=0.0D0
  DO ISD=1,NSD
    IA=IDXA(ISD)
    IB=IDXB(ISD)
    RHS%COEF(IA,IB)=X(ISD)
  END DO

  ! apply second series of operators exp(-Q_SI) exp(-Q_SA) exp(-Q_SI)

  ! exp(-Q_AI) = 1 - Sum_AI Q_AI E_AI
  SGM%COEF=0.0D0
  DO I=1,NI
    DO T=1,NA
      CALL WFN_EX1(-F%AI(T,I),NI+T,I,RHS,SGM)
    END DO
  END DO
  PSI1%COEF=SGM%COEF
  DO I=1,NI
    DO T=1,NA
      CALL WFN_EX1(-0.5D0*F%AI(T,I),NI+T,I,PSI1,SGM)
    END DO
  END DO
  RHS%COEF=RHS%COEF+SGM%COEF

  ! exp(-Q_SA) = 1 - Sum_SA Q_SA E_SA + ...
  SGM%COEF=0.0D0
  DO T=1,NA
    DO A=1,NS
      CALL WFN_EX1(-F%SA(A,T),NI+NA+A,NI+T,RHS,SGM)
    END DO
  END DO
  PSI1%COEF=SGM%COEF
  DO T=1,NA
    DO A=1,NS
      CALL WFN_EX1(-0.5D0*F%SA(A,T),NI+NA+A,NI+T,PSI1,SGM)
    END DO
  END DO
  RHS%COEF=RHS%COEF+SGM%COEF

  ! exp(-Q_SI) = 1 - Sum_SI Q_SI E_SI + ...
  SGM%COEF=0.0D0
  DO I=1,NI
    DO A=1,NS
      CALL WFN_EX1(-F%SI(A,I),NI+NA+A,I,RHS,SGM)
    END DO
  END DO
  PSI1%COEF=SGM%COEF
  DO I=1,NI
    DO A=1,NS
      CALL WFN_EX1(-0.5D0*F%SI(A,I),NI+NA+A,I,PSI1,SGM)
    END DO
  END DO
  RHS%COEF=RHS%COEF+SGM%COEF

  ! contract into X_i = <i|exp(-Q_IA)|RHS> and
  ! transform to adapt for the overcomplete basis
  DO ISD=1,NSD
    IA=IDXA(ISD)
    IB=IDXB(ISD)
    X(ISD)=RHS%COEF(IA,IB)
  END DO


  ! compute the energy

  ! method 1: E2 = <V0|X>
  E2 = DDOT_(NSD,V0,1,X,1)
  WRITE(*,'(1X,A32,F21.14)') 'E2 as <V0|X> = ', E2

  ! method 2: E2 = <0|Ĥ|Psi1>
  ! first, form the first-order wavefunction |Psi1>
  PSI1%COEF=0.0D0
  DO ISD=1,NSD
    IA=IDXA(ISD)
    IB=IDXB(ISD)
    PSI1%COEF(IA,IB)=X(ISD)
  END DO

  !WRITE(*,*) 'Final first-order wavefunction from UDL method:'
  !CALL WFN_PRINT(PSI1,1.0D-8)

  ! TAU still contains Ĥ|0>, so contract with |Psi1>
  E2 = DDOT_(PSI1%NDET,PSI1%COEF,1,TAU%COEF,1)
  WRITE(*,'(1X,A32,F21.14)') 'E2 as <Psi1|H|0> = ', E2

END SUBROUTINE CASPT2_UDL
