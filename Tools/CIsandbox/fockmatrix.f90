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
MODULE FOCKMATRIX
  ! implements a type for a blocked Fock matrix and provides functions for
  ! allocating/freeing, printing, and transforming the Fock matrix
  !
  ! Steven Vancoillie, November 2013, Lund
  USE ISO_FORTRAN_ENV
  IMPLICIT NONE

  TYPE FMAT
    REAL(REAL64), ALLOCATABLE :: &
       & II(:,:), IA(:,:), IS(:,:), &
       & AI(:,:), AA(:,:), AS(:,:), &
       & SI(:,:), SA(:,:), SS(:,:), &
       & DI(:), DA(:), DS(:)
  END type FMAT

CONTAINS

  SUBROUTINE FMAT_INIT(F,NI,NA,NS,FTOT,NMO)
    ! allocate Fock matrix blocks
    IMPLICIT NONE
    TYPE(FMAT), INTENT(OUT) :: F
    INTEGER, INTENT(IN) :: NI, NA, NS, NMO
    REAL(REAL64), INTENT(IN) :: FTOT(NMO,NMO)
    INTEGER :: NISTA, NIEND, NASTA, NAEND, NSSTA, NSEND
    ALLOCATE(F%II(NI,NI))
    ALLOCATE(F%IA(NI,NA))
    ALLOCATE(F%IS(NI,NS))
    ALLOCATE(F%AI(NA,NI))
    ALLOCATE(F%AA(NA,NA))
    ALLOCATE(F%AS(NA,NS))
    ALLOCATE(F%SI(NS,NI))
    ALLOCATE(F%SA(NS,NA))
    ALLOCATE(F%SS(NS,NS))
    ALLOCATE(F%DI(NI))
    ALLOCATE(F%DA(NA))
    ALLOCATE(F%DS(NS))
    NISTA=1
    NIEND=NI
    NASTA=NIEND+1
    NAEND=NIEND+NA
    NSSTA=NAEND+1
    NSEND=NAEND+NS
    F%II = FTOT(NISTA:NIEND,NISTA:NIEND)
    F%AI = FTOT(NASTA:NAEND,NISTA:NIEND)
    F%SI = FTOT(NSSTA:NSEND,NISTA:NIEND)
    F%IA = FTOT(NISTA:NIEND,NASTA:NAEND)
    F%AA = FTOT(NASTA:NAEND,NASTA:NAEND)
    F%SA = FTOT(NSSTA:NSEND,NASTA:NAEND)
    F%IS = FTOT(NISTA:NIEND,NSSTA:NSEND)
    F%AS = FTOT(NASTA:NAEND,NSSTA:NSEND)
    F%SS = FTOT(NSSTA:NSEND,NSSTA:NSEND)
  END SUBROUTINE FMAT_INIT

  SUBROUTINE FMAT_FREE(F)
    ! free Fock matrix blocks
    IMPLICIT NONE
    TYPE(FMAT) :: F
    DEALLOCATE(F%II)
    DEALLOCATE(F%IA)
    DEALLOCATE(F%IS)
    DEALLOCATE(F%AI)
    DEALLOCATE(F%AA)
    DEALLOCATE(F%AS)
    DEALLOCATE(F%SI)
    DEALLOCATE(F%SA)
    DEALLOCATE(F%SS)
    DEALLOCATE(F%DI)
    DEALLOCATE(F%DA)
    DEALLOCATE(F%DS)
  END SUBROUTINE FMAT_FREE

  SUBROUTINE FMAT_LOAD(F)
    USE ISO_FORTRAN_ENV
    IMPLICIT NONE
    TYPE(FMAT) :: F
    INTEGER :: I
    INTEGER :: IERR
    OPEN(UNIT=8,FILE='FOCK.DAT',STATUS='OLD',ACTION='READ',IOSTAT=IERR)
    IF (IERR.NE.0) STOP 'FMAT_LOAD: could not open Fock matrix file'
    DO I=1,SIZE(F%II,1)
      READ(8,*) F%II(I,:), F%IA(I,:), F%IS(I,:)
    END DO
    DO I=1,SIZE(F%AA,1)
      READ(8,*) F%AI(I,:), F%AA(I,:), F%AS(I,:)
    END DO
    DO I=1,SIZE(F%SS,1)
      READ(8,*) F%SI(I,:), F%SA(I,:), F%SS(I,:)
    END DO
    CLOSE(8)
  END SUBROUTINE FMAT_LOAD

  SUBROUTINE FMAT_COPY(F,FTOT,NMO)
    USE ISO_FORTRAN_ENV
    IMPLICIT NONE
    TYPE(FMAT), INTENT(IN) :: F
    INTEGER, INTENT(IN) :: NMO
    REAL(REAL64), INTENT(OUT) :: FTOT(NMO,NMO)
    INTEGER :: NI, NA, NS
    INTEGER :: NISTA, NIEND, NASTA, NAEND, NSSTA, NSEND
    NI=SIZE(F%II,1)
    NA=SIZE(F%AA,1)
    NS=SIZE(F%SS,1)
    NISTA=1
    NIEND=NI
    NASTA=NIEND+1
    NAEND=NIEND+NA
    NSSTA=NAEND+1
    NSEND=NAEND+NS
    FTOT(NISTA:NIEND,NISTA:NIEND) = F%II
    FTOT(NASTA:NAEND,NISTA:NIEND) = F%AI
    FTOT(NSSTA:NSEND,NISTA:NIEND) = F%SI
    FTOT(NISTA:NIEND,NASTA:NAEND) = F%IA
    FTOT(NASTA:NAEND,NASTA:NAEND) = F%AA
    FTOT(NSSTA:NSEND,NASTA:NAEND) = F%SA
    FTOT(NISTA:NIEND,NSSTA:NSEND) = F%IS
    FTOT(NASTA:NAEND,NSSTA:NSEND) = F%AS
    FTOT(NSSTA:NSEND,NSSTA:NSEND) = F%SS
  END SUBROUTINE FMAT_COPY

  SUBROUTINE FMAT_PRINT(F)
    USE ISO_FORTRAN_ENV
    IMPLICIT NONE
    TYPE(FMAT) :: F
    INTEGER :: I
    CHARACTER(LEN=64) :: FMT
    WRITE(*,*) 'FOCK MATRIX:'
    WRITE(FMT,*) '(F21.14,5X,2F21.14,2X,3F21.14,2X,2F21.14)'
    DO I=1,SIZE(F%II,1)
      WRITE(*,FMT) F%DI(I), F%II(I,:), F%IA(I,:), F%IS(I,:)
    END DO
    WRITE(*,*)
    DO I=1,SIZE(F%AA,1)
      WRITE(*,FMT) F%DA(I), F%AI(I,:), F%AA(I,:), F%AS(I,:)
    END DO
    WRITE(*,*)
    DO I=1,SIZE(F%SS,1)
      WRITE(*,FMT) F%DS(I), F%SI(I,:), F%SA(I,:), F%SS(I,:)
    END DO
    WRITE(*,*)
  END SUBROUTINE FMAT_PRINT

  SUBROUTINE FMAT_UDL(F)
    ! UDL decomposition of the Fock matrix.
    !
    ! The Fock matrix is brought into a diagonal form D, and upper and lower
    ! triangular matrices U and L. Intermediate changes in orbitals basis during
    ! the decomposition are collected in the orthogonal transformation matrices O,
    ! such that eventually, the original matrix F = O U D L O^T. The Fock matrix
    ! in the new basis is then F' = U D L, and F = O F' O^T. The storage of the
    ! new stuff overwrites the input blocked Fock matrix.
    !
    ! input: Fock matrix as blocks:
    !         FII FIA FIS
    !         FAI FAA FAS
    !         FSI FSA FSS
    ! output:
    !  - diagonal blocks FII, FAA, FSS contain the new orbital basis
    !  - new arrays DI, DS, DA contain the diagonal values that form the
    !    matrix D resulting from diagonalizing the diagonal blocks
    !  - off-diagonal blocks FIA, FIS, FAS contain the upper triangular
    !    blocks QIA, QIS, QAS that form the matrix U
    !  - off-diagonal blocks FAI, FSI, FSA contain the lower triangular
    !    blocks QAI, QSI, QSA that form the matrix L

    USE ISO_FORTRAN_ENV
    IMPLICIT NONE
    TYPE(FMAT), INTENT(INOUT) :: F

    INTEGER :: I, J
    INTEGER :: NI, NA, NS, NTOT
    INTEGER :: II, IA, IS
    INTEGER :: NISTA, NIEND, NASTA, NAEND, NSSTA, NSEND

    REAL(REAL64), ALLOCATABLE :: O(:,:), U(:,:), D(:,:), L(:,:)

    REAL(REAL64), ALLOCATABLE :: WORK(:)
    INTEGER :: NWORK, INFO

    NI=SIZE(F%DI)
    NA=SIZE(F%DA)
    NS=SIZE(F%DS)
    NTOT=NI+NA+NS

    ! allocate temporary scratch
    NWORK=1+MAX(MAX(NI,NA),NS)**2
    ALLOCATE(WORK(NWORK))

    ! step 1: diagonalize the secondary-secondary subblock
    call dsyev_('V','U',NS,F%SS,NS,F%DS,WORK,NWORK,INFO)
    IF (INFO.NE.0) STOP 'Error: diagonalization of FSS failed'

    ! adjust off-diagonal blocks to basis change
    F%IS = MATMUL(F%IS,F%SS)
    F%AS = MATMUL(F%AS,F%SS)
    F%SI = TRANSPOSE(F%IS)
    F%SA = TRANSPOSE(F%AS)

    ! construct new QIS/QSI blocks
    DO IS=1,NS
      F%IS(:,IS) = F%IS(:,IS)/F%DS(IS)
    END DO
    F%II = F%II - MATMUL(F%IS,F%SI)
    F%SI = TRANSPOSE(F%IS)

    ! adjust off-diagonal blocks to elimination
    F%AI = F%AI - MATMUL(F%AS,F%SI)
    F%IA = TRANSPOSE(F%AI)

    ! construct new QAS/QSA blocks
    DO IS=1,NS
      F%AS(:,IS) = F%AS(:,IS)/F%DS(IS)
    END DO
    F%AA = F%AA - MATMUL(F%AS,F%SA)
    F%SA = TRANSPOSE(F%AS)

    ! step 2: diagonalize the active-active subblock
    call dsyev_('V','U',NA,F%AA,NA,F%DA,WORK,NWORK,INFO)
    IF (INFO.NE.0) STOP 'Error: diagonalization of FAA failed'

    ! adjust off-diagonal blocks to basis change
    F%IA = MATMUL(F%IA,F%AA)
    F%SA = MATMUL(F%SA,F%AA)
    F%AI = TRANSPOSE(F%IA)
    F%AS = TRANSPOSE(F%SA)

    ! construct new QIA/QAI blocks
    DO IA=1,NA
      F%IA(:,IA) = F%IA(:,IA)/F%DA(IA)
    END DO
    F%II = F%II - MATMUL(F%IA,F%AI)
    F%AI = TRANSPOSE(F%IA)

    ! step 3: diagonalize the inactive-inactive subblock
    call dsyev_('V','U',NI,F%II,NI,F%DI,WORK,NWORK,INFO)
    IF (INFO.NE.0) STOP 'Error: diagonalization of FII failed'

    ! adjust off-diagonal blocks to basis change
    F%AI = MATMUL(F%AI,F%II)
    F%SI = MATMUL(F%SI,F%II)
    F%IA = TRANSPOSE(F%AI)
    F%IS = TRANSPOSE(F%SI)

#ifdef _DEBUG_
    ! sanity test: we should be able to reconstruct the original fock matrix!!
    ALLOCATE(O(NTOT,NTOT))
    ALLOCATE(U(NTOT,NTOT))
    ALLOCATE(D(NTOT,NTOT))
    ALLOCATE(L(NTOT,NTOT))

    O = 0.0D0
    U = 0.0D0
    D = 0.0D0
    L = 0.0D0
    DO I=1,NTOT
      U(I,I) = 1.0D0
      L(I,I) = 1.0D0
    END DO

    NISTA=1
    NIEND=NI
    NASTA=NIEND+1
    NAEND=NIEND+NA
    NSSTA=NAEND+1
    NSEND=NAEND+NS

    O(NISTA:NIEND,NISTA:NIEND) = F%II
    O(NASTA:NAEND,NASTA:NAEND) = F%AA
    O(NSSTA:NSEND,NSSTA:NSEND) = F%SS

    WRITE(*,*) 'Orthogonal basis transform:'
    DO I=1,NTOT
      WRITE(*,'(10F10.6)') O(I,:)
    END DO

    U(NISTA:NIEND,NASTA:NAEND) = F%IA
    U(NISTA:NIEND,NSSTA:NSEND) = F%IS
    U(NASTA:NAEND,NSSTA:NSEND) = F%AS

    WRITE(*,*) 'Upper:'
    DO I=1,NTOT
      WRITE(*,'(10F10.6)') U(I,:)
    END DO

    L(NASTA:NAEND,NISTA:NIEND) = F%AI
    L(NSSTA:NSEND,NISTA:NIEND) = F%SI
    L(NSSTA:NSEND,NASTA:NAEND) = F%SA

    WRITE(*,*) 'Lower:'
    DO I=1,NTOT
      WRITE(*,'(10F10.6)') L(I,:)
    END DO

    DO II=1,NI
      D(II,II) = F%DI(II)
    END DO
    DO IA=1,NA
      D(NASTA+IA-1,NASTA+IA-1) = F%DA(IA)
    END DO
    DO IS=1,NS
      D(NSSTA+IS-1,NSSTA+IS-1) = F%DS(IS)
    END DO

    WRITE(*,*) 'Diagonal:'
    DO I=1,NTOT
      WRITE(*,'(10F10.6)') D(I,:)
    END DO

    O = MATMUL(MATMUL(MATMUL(MATMUL(O,U),D),L),TRANSPOSE(O))

    WRITE(*,*) 'Back to Fock matrix???'
    DO I=1,NTOT
      WRITE(*,'(10F10.6)') O(I,:)
    END DO
    DEALLOCATE(O,U,D,L)
#endif

  END SUBROUTINE FMAT_UDL

END MODULE
