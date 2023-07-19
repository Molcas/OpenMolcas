************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2023, Ignacio Fdez. Galvan                             *
************************************************************************

      ! Sort Dyson orbitals according to symmetry so that they can be
      ! written with WrVec

      SUBROUTINE WRVEC_DYSON(filename,LUNIT,NSYM,NBAS,ORBNUM,CMO,AMPS,
     &                       DYSEN,TITLE,NZ)

      USE stdalloc, ONLY: mma_allocate, mma_deallocate
      USE Constants, ONLY: Zero, One

      IMPLICIT NONE
      CHARACTER(LEN=*) :: filename, TITLE
      INTEGER :: LUNIT, NSYM, NBAS(NSYM), ORBNUM, NZ
      REAL*8 :: CMO(NZ,ORBNUM), AMPS(ORBNUM), DYSEN(ORBNUM)
      INTEGER :: DUMMY(7,8), I, J, NB(0:NSYM), NBAS_(NSYM), NBT,
     &           NORB(NSYM), NSYM_, OFF, ORBSYM(ORBNUM)
      REAL*8 :: RSUM
      LOGICAL :: DODESYM
      REAL*8, ALLOCATABLE :: DESYM(:,:), REORD(:)

      ! First count how many orbitals in each symmetry
      NB(0) = 1
      DO I=1,NSYM
        NB(I) = NB(I-1)+NBAS(I)
      END DO
      DODESYM = .FALSE.
      NORB(:) = 0
      ORBSYM(:) = 0
      outer: DO I=1,ORBNUM
        DO J=1,NSYM
          RSUM = SUM(ABS(CMO(NB(J-1):NB(J)-1,I)))
          IF (RSUM > Zero) THEN
            ! If there is any orbital with mixture of symmetries,
            ! we have to desymmetrize the whole thing
            IF (ORBSYM(I) /= 0) THEN
              DODESYM = .TRUE.
              EXIT outer
            END IF
            ORBSYM(I) = J
            NORB(J) = NORB(J)+1
          END IF
        END DO
      END DO outer
      IF (DODESYM) THEN
        NSYM_ = 1
        NBAS_(1) = SUM(NBAS(1:NSYM))
        NB(1) = NBAS_(1)
        NORB(1) = ORBNUM
        ORBSYM(:) = 1
      ELSE
        NSYM_ = NSYM
        NBAS_(:) = NBAS(:)
      END IF
      NBT = 0
      DO I=1,NSYM_
        NBT = NBT+NORB(I)*NBAS_(I)
      END DO
      CALL mma_allocate(REORD,NBT,Label='REORD')
      IF (DODESYM) THEN
        ! Here do a plain desymmetrization
        NBT = NBAS_(1)
        CALL mma_allocate(DESYM,NBT,NBT,Label='DESYM')
        CALL get_dArray('SM',DESYM,NBT**2)
        CALL DGEMM_('N','N',NBT,NORB(1),NBT,One,DESYM,NBT,CMO,NBT,Zero,
     &              REORD,NBT)
        CALL mma_deallocate(DESYM)
      ELSE
        ! Here distribute the coefficients
        OFF = 0
        DO I=1,NSYM_
          DO J=1,ORBNUM
            IF (ORBSYM(J) /= I) CYCLE
            REORD(OFF+1:OFF+NBAS_(I)) = CMO(NB(I-1):NB(I)-1,J)
            OFF = OFF+NBAS_(I)
          END DO
        END DO
      END IF

      ! And call WrVec with the reordered data
      CALL WRVEC(filename,LUNIT,'COE',NSYM_,NBAS_,NORB,REORD,AMPS,
     &           DYSEN,DUMMY,TITLE)

      CALL mma_deallocate(REORD)

      RETURN

      END
