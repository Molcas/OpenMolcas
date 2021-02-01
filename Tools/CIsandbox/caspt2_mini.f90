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
PROGRAM CASPT2_MINI
  ! A mini-caspt2 program to test the new UDL decompostition in a flexible
  ! environment.  This initial version is based on a Molcas 7.9 calculation of
  ! the 'crazy' system BeHe2Li (without IPEA shift!).

  ! Steven Vancoillie, November-December 2013, Lund

  USE ISO_FORTRAN_ENV, ONLY: REAL64
  USE SECOND_QUANTIZATION
  USE WAVEFUNCTION
  USE DENSITY
  USE FOCKMATRIX
  USE ORBINT
  IMPLICIT NONE

  ! orbital indices
  INTEGER :: P, Q, R, S

  ! the wavefunction
  TYPE(WFN) :: PSI
  INTEGER :: NEL, NORB, MULT

  ! Fock matrix
  TYPE(FMAT) :: F
  INTEGER :: NI, NA, NS, NAO, NMO
  REAL(REAL64), ALLOCATABLE :: FTOT(:,:)

  ! orbint stuff
  REAL(REAL64), ALLOCATABLE :: CMO(:,:)
  REAL(REAL64), ALLOCATABLE :: ONEINT(:,:)
  REAL(REAL64), ALLOCATABLE :: TWOINT(:,:,:,:)

  ! the one- and two- electron density matrices
  REAL(REAL64), ALLOCATABLE :: D1(:,:)
  REAL(REAL64), ALLOCATABLE :: D2(:,:,:,:)

  ! extra info
  INTEGER :: NTOTEL

  ! core energy
  REAL(REAL64) :: ECORE

  ! initialize linalg library
  CALL INIT_LINALG

  ! initialize binomial tables needed for ranking
  CALL RANK_INIT

  !################################################################
  !#### Set up the wavefunction
  !################################################################

  ! initialize values to match the 'crazy' BeHeLi system
  NEL = 3
  NORB = 3
  MULT = 2

  NI=2
  NA=3
  NS=2
  NMO=NI+NA+NS
  NAO=38

  ECORE=-16.704241392292605D0

  ! initialize the wavefunction
  CALL WFN_INIT(PSI,NEL,NORB,MULT)

  ! set coefficients to the rasscf coefficients
  ! of the LiHe2Be calculations
  PSI%COEF    =  0.0D0
  PSI%COEF(1,1) =  1.d0 * ( 0.984615412797048d0) * (-1.d0)
  PSI%COEF(1,2) = -1.d0 * ( 0.116525050204646d-1)
  PSI%COEF(2,1) =  1.d0 * (-0.332209883242451d-12) * (-1.d0)
  PSI%COEF(2,2) =  SQRT(1.d0/2.d0) * (-0.332387701248167d-11) * (-1.d0)
  PSI%COEF(3,1) = -SQRT(1.d0/2.d0) * (-0.332387701248167d-11)
  PSI%COEF(3,2) =  1.d0 * ( 0.250355292052973d-12) * (-1.d0)
  PSI%COEF(1,3) = -SQRT(2.d0/3.d0) * (-0.194099496829364d-11)
  PSI%COEF(2,2) = PSI%COEF(2,2) + SQRT(1.d0/6.d0) * (-0.194099496829364d-11) *(-1.d0)
  PSI%COEF(3,1) = PSI%COEF(3,1) + SQRT(1.d0/6.d0) * (-0.194099496829364d-11)
  PSI%COEF(2,3) = -1.d0 * (-0.723233381792119d-1)
  PSI%COEF(3,3) = -1.d0 * ( 0.158638087368269d0)

  ! normalize the CI coefficients if needed
  !CALL WFN_NORMALIZE(PSI)

  ! print the wavefunction, takes CI coefficient treshold
  CALL WFN_PRINT(PSI,0.0D0)

  !################################################################
  !#### Compute density matrices
  !################################################################

  ALLOCATE(D1(NORB,NORB))
  CALL D1_ANN(PSI,D1)

  ALLOCATE(D2(NORB,NORB,NORB,NORB))
  CALL D2_ANN(PSI,D2)

  !################################################################
  !#### Set up the MO integrals
  !################################################################

  ALLOCATE(CMO(NAO,NMO))
  ALLOCATE(ONEINT(NMO,NMO))
  ALLOCATE(TWOINT(NMO,NMO,NMO,NMO))

  CALL ORBINT_LOAD(CMO,ONEINT,TWOINT)

  !################################################################
  !#### Compute energy
  !################################################################

  CALL WFN_ENERGY(NI,NA,NS,D1,D2,ONEINT,TWOINT)

  !################################################################
  !#### Set up the Fock matrix
  !################################################################

  ! compute Fock matrix
  ALLOCATE(FTOT(NMO,NMO))
  DO P=1,NMO
    DO Q=1,NMO
      FTOT(P,Q)=ONEINT(P,Q)
      DO R=1,NI
        FTOT(P,Q)=FTOT(P,Q)+2.0D0*(TWOINT(P,Q,R,R)-0.5D0*TWOINT(P,R,Q,R))
      END DO
      DO R=NI+1,NI+NA
        DO S=NI+1,NI+NA
          FTOT(P,Q)=FTOT(P,Q)+D1(R-NI,S-NI)*(TWOINT(P,Q,R,S)-0.5D0*TWOINT(P,R,Q,S))
        END DO
      END DO
    END DO
  END DO

  CALL FMAT_INIT(F,NI,NA,NS,FTOT,NMO)
  F%DI=0.0D0
  F%DA=0.0D0
  F%DS=0.0D0

  CALL FMAT_PRINT(F)

  !################################################################
  !#### Regular CASPT2
  !################################################################

  ! Solve CASPT2 by solving a linear system of equations
  !WRITE(*,*) 'CASPT2 Ax=b'
  CALL CASPT2_AXB(PSI,F,NMO,ONEINT,TWOINT)

  !################################################################
  !#### Direct CASPT2
  !################################################################

  WRITE(*,*) 'CASPT2 UDL'
  CALL CASPT2_UDL(PSI,F,NAO,NMO,CMO,ONEINT,TWOINT)

  !################################################################
  !#### Full CI
  !################################################################

  NTOTEL=2*NI+NEL
  !CALL FULLCI(NTOTEL,NMO,MULT,ONEINT,TWOINT,ECORE)

END PROGRAM CASPT2_MINI
