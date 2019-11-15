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
SUBROUTINE FULLCI(NEL,NORB,MULT,ONEINT,TWOINT,ECORE)
  ! Computes the CASPT2 second-order energy of the wavefunction PSI, provided
  ! that the full Fock matrix and necessary 1- and 2-el integrals are given.
  USE ISO_FORTRAN_ENV, ONLY: REAL64
  USE SECOND_QUANTIZATION
  USE WAVEFUNCTION
  USE DENSITY
  USE FOCKMATRIX
  USE ORBINT
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: NEL, NORB, MULT
  REAL(REAL64), INTENT(IN) :: ONEINT(NORB,NORB), TWOINT(NORB,NORB,NORB,NORB)
  REAL(REAL64), INTENT(IN) :: ECORE

  TYPE(WFN) :: PSI, SGM
  REAL(REAL64), ALLOCATABLE :: D1(:,:), D2(:,:,:,:)

  REAL(REAL64), ALLOCATABLE :: H(:,:,:,:), E(:)

  REAL(REAL64), ALLOCATABLE :: WORK(:)
  INTEGER, ALLOCATABLE :: IWORK(:)
  INTEGER :: LWORK, LIWORK, ISUPPZ(2), INFO
  REAL(REAL64) :: VL, VU
  INTEGER :: IL, IU

  INTEGER :: IA, IB, DETA, DETB, NDET
  INTEGER :: P, Q, R, S

  ! initialize the wavefunction

  CALL WFN_INIT(PSI,NEL,NORB,MULT)
  CALL WFN_INIT(SGM,NEL,NORB,MULT)

  ! construct the Full-CI hamiltonian

  ALLOCATE(H(PSI%NDETA,PSI%NDETB,PSI%NDETA,PSI%NDETB))

  DETB=LEX_INIT(PSI%NELB,PSI%NORB)
  DO IB=1,PSI%NDETB
    DETA=LEX_INIT(PSI%NELA,PSI%NORB)
    DO IA=1,PSI%NDETA

      SGM%COEF=0.0D0
      ! operate with sum over h_pq E_pq
      DO P=1,NORB
        DO Q=1,NORB
          CALL DET_EX1(ONEINT(P,Q),P,Q,DETA,DETB,SGM)
          DO R=1,NORB
            DO S=1,NORB
              CALL DET_EX2(0.5D0*TWOINT(P,Q,R,S),P,Q,R,S,DETA,DETB,SGM)
            END DO
          END DO
        END DO
      END DO

      H(:,:,IA,IB)=SGM%COEF
      H(IA,IB,IA,IB)=H(IA,IB,IA,IB)+ECORE

      DETA=LEX_NEXT(DETA)
    END DO
    DETB=LEX_NEXT(DETB)
  END DO

  ! diagonalize the entire hamiltonian
  NDET=PSI%NDETA*PSI%NDETB
  LWORK=NDET**2
  LIWORK=10*NDET
  ALLOCATE(WORK(LWORK))
  ALLOCATE(IWORK(LIWORK))
  ALLOCATE(E(NDET))
  IL=1
  IU=1

  ! find all eigenvalues/eigenvectors of H
  !call dsyev_('V','U',NCI,H,NCI,SPEC,WORK,NWORK,INFO)

  ! find lowest eigenvalue/eigenvector of H
  call dsyevr_('V','I','U',NDET,H,NDET,VL,VU,IL,IU,0.0D0,1,E,PSI%COEF,NDET,ISUPPZ,WORK,LWORK,IWORK,LIWORK,INFO)
  IF (INFO.NE.0) STOP 'FULLCI: diagonalization of FCI hamiltonian failed'

  WRITE(*,'(1X,A,F21.14)') 'Full-CI GS energy:', E(1)

  CALL WFN_PRINT(PSI,0.05D0)

  ALLOCATE(D1(NORB,NORB))
  ALLOCATE(D2(NORB,NORB,NORB,NORB))
  CALL D1_ANN(PSI,D1)
  CALL D2_ANN(PSI,D2)
  
  CALL WFN_ENERGY(0,NORB,0,D1,D2,ONEINT,TWOINT)

END SUBROUTINE FULLCI
