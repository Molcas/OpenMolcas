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
      SUBROUTINE SPINORB(D,CMO,OCC,kroot)
      use stdalloc, only: mma_allocate, mma_deallocate
      use PrintLevel, only: DEBUG
      use output_ras, only: LF,IPRLOC
      use general_data, only: NSYM,NASH,NBAS,NFRO,NISH

!
!     Purpose: diagonalize the spin density matrix (D) to
!     obtain the eigenvectors (EVEC) and the eigenvalues (EVAL).
!     Then the natural spinorbitals (CMONSO) are computed
!     (only active).
!
      IMPLICIT None
      Real*8 D(*),CMO(*),OCC(*)
      Integer :: KROOT

      Character(LEN=16):: ROUTINE='SPINORB '
      Real*8, Allocatable:: W1(:,:), W2(:,:)
      Integer :: I, IPCMO, IPDEN, IPOCC, iPrLev, iSym, NA, NB, NF, NI
      Integer :: IDIAG
!
!
! Local print level (if any)
      IPRLEV=IPRLOC(6)
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF

      IPDEN=1
      IPCMO=1
      IPOCC=1
      DO ISYM=1,NSYM
        NB=NBAS(ISYM)
        NF=NFRO(ISYM)
        NI=NISH(ISYM)
        IF ( NB.NE.0 ) THEN
          NA=NASH(ISYM)
          IF ( NA.NE.0 ) THEN
            CALL mma_allocate(W1,NA,NA,Label='W1')
            CALL mma_allocate(W2,NB,NA,Label='W2')
            W1(:,:)=0.0D0
            CALL DCOPY_(NA,[1.0D0],0,W1,NA+1)
            CALL Jacob(D(IPDEN),W1,NA,NA)
            IDIAG=0
            DO I=1,NA
              IDIAG=IDIAG+I
              OCC(IPOCC+NF+NI+I-1)=D(IPDEN+IDIAG-1)
            END DO
            CALL DGEMM_('N','N',                                        &
     &                  NB,NA,NA,                                       &
     &                  1.0d0,CMO(IPCMO+(NF+NI)*NB),NB,                 &
     &                        W1,NA,                                    &
     &                  0.0d0,W2,NB)
            CALL DCOPY_(NA*NB,W2,1,CMO(IPCMO+(NF+NI)*NB),1)
            Call mma_deallocate(W2)
            Call mma_deallocate(W1)
            IPDEN=IPDEN+NA*(NA+1)/2
          END IF
          IPCMO=IPCMO+NB*NB
          IPOCC=IPOCC+NB
        END IF
      END DO
!
!
! Avoid unused argument warnings
#ifdef _WARNING_WORKAROUND_
      IF (.FALSE.) CALL Unused_integer(kroot)
#endif
      END SUBROUTINE SPINORB
