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

  USE ISO_FORTRAN_ENV, ONLY: REAL64
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
    REAL(REAL64), ALLOCATABLE :: TMP1(:,:), TMP2(:,:)
    INTEGER :: NI, NA, NS, NMO
    INTEGER :: P, Q

    NI=SIZE(OII,1)
    NA=SIZE(OAA,1)
    NS=SIZE(OSS,1)

    NMO=SIZE(CMO,2)

    ! orbitals
    ALLOCATE(TMP1(NMO,NI))
    TMP1(:,:) = MATMUL(CMO(:,1:NI),OII)
    CMO(:,1:NI) = TMP1
    DEALLOCATE(TMP1)
    ALLOCATE(TMP1(NMO,NA))
    TMP1(:,:) = MATMUL(CMO(:,NI+1:NI+NA),OAA)
    CMO(:,NI+1:NI+NA) = TMP1
    DEALLOCATE(TMP1)
    ALLOCATE(TMP1(NMO,NS))
    TMP1(:,:) = MATMUL(CMO(:,NI+NA+1:NI+NA+NS),OSS)
    CMO(:,NI+NA+1:NI+NA+NS) = TMP1
    DEALLOCATE(TMP1)

    ! 1-el integral transforms
    ! inactive
    ALLOCATE(TMP1(NMO,NI),TMP2(NI,NMO))
    TMP1(:,:) = MATMUL(ONEINT(:,1:NI),OII)
    ONEINT(:,1:NI) = TMP1
    TMP2(:,:) = MATMUL(TRANSPOSE(OII),ONEINT(1:NI,:))
    ONEINT(1:NI,:) = TMP2
    DEALLOCATE(TMP1,TMP2)
    ! active
    ALLOCATE(TMP1(NMO,NA),TMP2(NA,NMO))
    TMP1(:,:) = MATMUL(ONEINT(:,NI+1:NI+NA),OAA)
    ONEINT(:,NI+1:NI+NA) = TMP1
    TMP2(:,:) = MATMUL(TRANSPOSE(OAA),ONEINT(NI+1:NI+NA,:))
    ONEINT(NI+1:NI+NA,:) = TMP2
    DEALLOCATE(TMP1,TMP2)
    ! virtual
    ALLOCATE(TMP1(NMO,NS),TMP2(NS,NMO))
    TMP1(:,:) = MATMUL(ONEINT(:,NI+NA+1:NI+NA+NS),OSS)
    ONEINT(:,NI+NA+1:NI+NA+NS) = TMP1
    TMP2(:,:) = MATMUL(TRANSPOSE(OSS),ONEINT(NI+NA+1:NI+NA+NS,:))
    ONEINT(NI+NA+1:NI+NA+NS,:) = TMP2
    DEALLOCATE(TMP1,TMP2)

    ! 2-el integral transforms
    ! inactive
    ALLOCATE(TMP1(NMO,NI),TMP2(NI,NMO))
    DO P=1,NMO
      DO Q=1,NMO
        TMP1(:,:) = MATMUL(TWOINT(:,P,Q,1:NI),OII)
        TWOINT(:,P,Q,1:NI) = TMP1
        TMP2(:,:) = MATMUL(TRANSPOSE(OII),TWOINT(1:NI,P,Q,:))
        TWOINT(1:NI,P,Q,:) = TMP2
      END DO
    END DO
    DO P=1,NMO
      DO Q=1,NMO
        TMP1(:,:) = MATMUL(TWOINT(P,:,1:NI,Q),OII)
        TWOINT(P,:,1:NI,Q) = TMP1
        TMP2(:,:) = MATMUL(TRANSPOSE(OII),TWOINT(P,1:NI,:,Q))
        TWOINT(P,1:NI,:,Q) = TMP2
      END DO
    END DO
    DEALLOCATE(TMP1,TMP2)
    ! active
    ALLOCATE(TMP1(NMO,NA),TMP2(NA,NMO))
    DO P=1,NMO
      DO Q=1,NMO
        TMP1(:,:) = MATMUL(TWOINT(:,P,Q,NI+1:NI+NA),OAA)
        TWOINT(:,P,Q,NI+1:NI+NA) = TMP1
        TMP2(:,:) = MATMUL(TRANSPOSE(OAA),TWOINT(NI+1:NI+NA,P,Q,:))
        TWOINT(NI+1:NI+NA,P,Q,:) = TMP2
      END DO
    END DO
    DO P=1,NMO
      DO Q=1,NMO
        TMP1(:,:) = MATMUL(TWOINT(P,:,NI+1:NI+NA,Q),OAA)
        TWOINT(P,:,NI+1:NI+NA,Q) = TMP1
        TMP2(:,:) = MATMUL(TRANSPOSE(OAA),TWOINT(P,NI+1:NI+NA,:,Q))
        TWOINT(P,NI+1:NI+NA,:,Q) = TMP2
      END DO
    END DO
    DEALLOCATE(TMP1,TMP2)
    ! virtual
    ALLOCATE(TMP1(NMO,NS),TMP2(NS,NMO))
    DO P=1,NMO
      DO Q=1,NMO
        TMP1(:,:) = MATMUL(TWOINT(:,P,Q,NI+NA+1:NI+NA+NS),OSS)
        TWOINT(:,P,Q,NI+NA+1:NI+NA+NS) = TMP1
        TMP2(:,:) = MATMUL(TRANSPOSE(OSS),TWOINT(NI+NA+1:NI+NA+NS,P,Q,:))
        TWOINT(NI+NA+1:NI+NA+NS,P,Q,:) = TMP2
      END DO
    END DO
    DO P=1,NMO
      DO Q=1,NMO
        TMP1(:,:) = MATMUL(TWOINT(P,:,NI+NA+1:NI+NA+NS,Q),OSS)
        TWOINT(P,:,NI+NA+1:NI+NA+NS,Q) = TMP1
        TMP2(:,:) = MATMUL(TRANSPOSE(OSS),TWOINT(P,NI+NA+1:NI+NA+NS,:,Q))
        TWOINT(P,NI+NA+1:NI+NA+NS,:,Q) = TMP2
      END DO
    END DO
    DEALLOCATE(TMP1,TMP2)
  END SUBROUTINE ORBINT_TRANSFORM

END MODULE ORBINT
