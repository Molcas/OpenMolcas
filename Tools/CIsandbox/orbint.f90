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
MODULE ORBINT
  ! Implements a way to store and transform the molecular orbital coefficients
  ! and the 1-and 2-electron MO integrals.
  !
  ! Steven Vancoillie, November 2013, Lund

  USE ISO_FORTRAN_ENV
  IMPLICIT NONE

CONTAINS

  SUBROUTINE ORBINT_LOAD(CMO,ONEINT,TWOINT)
    IMPLICIT NONE
    REAL(REAL64), INTENT(OUT) :: CMO(:,:), ONEINT(:,:), TWOINT(:,:,:,:)
    INTEGER :: NAO, NMO
    INTEGER :: I, J, P
    INTEGER :: IERR

    NAO=SIZE(CMO,1)
    NMO=SIZE(CMO,2)
    
    ! MO coefficients
    OPEN(UNIT=8,FILE='MOS.DAT',STATUS='OLD',ACTION='READ',IOSTAT=IERR)
    IF (IERR.NE.0) STOP 'ORBINT_INIT: could not open MO coefficients file'
    DO I=1,NAO
      READ(8,*) CMO(I,:)
    END DO
    CLOSE(8)

    ! element ONEINT(P,Q) stores h_pq
    OPEN(UNIT=8,FILE='ONEINT.DAT',STATUS='OLD',ACTION='READ',IOSTAT=IERR)
    IF (IERR.NE.0) STOP 'ORBINT_INIT: could not open 1-el integral file'
    DO I=1,NMO
      READ(8,*) ONEINT(I,1:I)
      ONEINT(1:I,I)=ONEINT(I,1:I)
    END DO
    CLOSE(8)

    ! element TWOINT(P,Q,I,J) stores g_pqij (in file as ij|pq)
    OPEN(UNIT=8,FILE='TWOINT.DAT',STATUS='OLD',ACTION='READ',IOSTAT=IERR)
    IF (IERR.NE.0) STOP 'ORBINT_INIT: could not open 1-el integral file'
    DO I=1,NMO
      DO J=1,NMO
        DO P=1,NMO
          READ(8,*) TWOINT(P,:,I,J)
        END DO
      END DO
    END DO
    CLOSE(8)
  END SUBROUTINE ORBINT_LOAD
  
  SUBROUTINE ORBINT_TRANSFORM(CMO,ONEINT,TWOINT,OII,OAA,OSS)
    IMPLICIT NONE
    REAL(REAL64), INTENT(IN) :: OII(:,:), OAA(:,:), OSS(:,:)
    REAL(REAL64), INTENT(INOUT) :: CMO(:,:), ONEINT(:,:), TWOINT(:,:,:,:)
    INTEGER :: NI, NA, NS, NMO
    INTEGER :: P, Q

    NI=SIZE(OII,1)
    NA=SIZE(OAA,1)
    NS=SIZE(OSS,1)

    NMO=SIZE(CMO,2)

    ! orbitals
    CMO(:,1:NI) = MATMUL(CMO(:,1:NI),OII)
    CMO(:,NI+1:NI+NA) = MATMUL(CMO(:,NI+1:NI+NA),OAA)
    CMO(:,NI+NA+1:NI+NA+NS) = MATMUL(CMO(:,NI+NA+1:NI+NA+NS),OSS)

    ! 1-el integral transforms
    ! inactive
    ONEINT(:,1:NI) = &
         & MATMUL(ONEINT(:,1:NI),OII)
    ONEINT(1:NI,:) = &
         & MATMUL(TRANSPOSE(OII),ONEINT(1:NI,:))
    ! active
    ONEINT(:,NI+1:NI+NA) = &
         & MATMUL(ONEINT(:,NI+1:NI+NA),OAA)
    ONEINT(NI+1:NI+NA,:) = &
         & MATMUL(TRANSPOSE(OAA),ONEINT(NI+1:NI+NA,:))
    ! virtual
    ONEINT(:,NI+NA+1:NI+NA+NS) = &
         & MATMUL(ONEINT(:,NI+NA+1:NI+NA+NS),OSS)
    ONEINT(NI+NA+1:NI+NA+NS,:) = &
         & MATMUL(TRANSPOSE(OSS),ONEINT(NI+NA+1:NI+NA+NS,:))

    ! 2-el integral transforms
    ! inactive
    DO P=1,NMO
      DO Q=1,NMO
        TWOINT(:,P,Q,1:NI) = &
             & MATMUL(TWOINT(:,P,Q,1:NI),OII)
        TWOINT(1:NI,P,Q,:) = &
             & MATMUL(TRANSPOSE(OII),TWOINT(1:NI,P,Q,:))
      END DO
    END DO
    DO P=1,NMO
      DO Q=1,NMO
        TWOINT(P,:,1:NI,Q) = &
             & MATMUL(TWOINT(P,:,1:NI,Q),OII)
        TWOINT(P,1:NI,:,Q) = &
             & MATMUL(TRANSPOSE(OII),TWOINT(P,1:NI,:,Q))
      END DO
    END DO
    ! active
    DO P=1,NMO
      DO Q=1,NMO
        TWOINT(:,P,Q,NI+1:NI+NA) = &
             & MATMUL(TWOINT(:,P,Q,NI+1:NI+NA),OAA)
        TWOINT(NI+1:NI+NA,P,Q,:) = &
             & MATMUL(TRANSPOSE(OAA),TWOINT(NI+1:NI+NA,P,Q,:))
      END DO
    END DO
    DO P=1,NMO
      DO Q=1,NMO
        TWOINT(P,:,NI+1:NI+NA,Q) = &
             & MATMUL(TWOINT(P,:,NI+1:NI+NA,Q),OAA)
        TWOINT(P,NI+1:NI+NA,:,Q) = &
             & MATMUL(TRANSPOSE(OAA),TWOINT(P,NI+1:NI+NA,:,Q))
      END DO
    END DO
    ! virtual
    DO P=1,NMO
      DO Q=1,NMO
        TWOINT(:,P,Q,NI+NA+1:NI+NA+NS) = &
             & MATMUL(TWOINT(:,P,Q,NI+NA+1:NI+NA+NS),OSS)
        TWOINT(NI+NA+1:NI+NA+NS,P,Q,:) = &
             & MATMUL(TRANSPOSE(OSS),TWOINT(NI+NA+1:NI+NA+NS,P,Q,:))
      END DO
    END DO
    DO P=1,NMO
      DO Q=1,NMO
        TWOINT(P,:,NI+NA+1:NI+NA+NS,Q) = &
             & MATMUL(TWOINT(P,:,NI+NA+1:NI+NA+NS,Q),OSS)
        TWOINT(P,NI+NA+1:NI+NA+NS,:,Q) = &
             & MATMUL(TRANSPOSE(OSS),TWOINT(P,NI+NA+1:NI+NA+NS,:,Q))
      END DO
    END DO
  END SUBROUTINE ORBINT_TRANSFORM

END MODULE ORBINT
