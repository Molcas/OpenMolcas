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
      USE Constants, ONLY: Zero

      IMPLICIT NONE
      CHARACTER(LEN=*) :: filename, TITLE
      INTEGER :: LUNIT, NSYM, NBAS(NSYM), ORBNUM, NZ
      REAL*8 :: CMO(NZ,ORBNUM), AMPS(ORBNUM), DYSEN(ORBNUM)
      INTEGER :: I, J, NB(0:NSYM), NBT, NORB(NSYM), OFF, ORBSYM(ORBNUM)
      REAL*8 :: DUMMY(1), RSUM
      REAL*8, ALLOCATABLE :: REORD(:)
      REAL*8, EXTERNAL :: dDot_

      ! First count how many orbitals in each symmetry
      NB(0) = 1
      DO I=1,NSYM
        NB(I) = NB(I-1)+NBAS(I)
      END DO
      NORB(:) = 0
      ORBSYM(:) = 0
      DO I=1,ORBNUM
        DO J=1,NSYM
          RSUM = SUM(ABS(CMO(NB(J-1):NB(J)-1,I)))
          IF (RSUM > Zero) THEN
            IF (ORBSYM(I) /= 0) THEN
              WRITE(6,100) TRIM(filename)
              RETURN
            END IF
            ORBSYM(I) = J
            NORB(J) = NORB(J)+1
          END IF
        END DO
      END DO
      ! Now distribute the coefficients
      NBT = 0
      DO I=1,NSYM
        NBT = NBT+NORB(I)*NBAS(I)
      END DO
      CALL mma_allocate(REORD,NBT,Label='REORD')
      OFF = 0
      DO I=1,NSYM
        DO J=1,ORBNUM
          IF (ORBSYM(J) /= I) CYCLE
          REORD(OFF+1:OFF+NBAS(I)) = CMO(NB(I-1):NB(I)-1,J)
          OFF = OFF+NBAS(I)
        END DO
      END DO

      ! And call WrVec with the reordered data
      CALL WRVEC(filename,LUNIT,'COE',NSYM,NBAS,NORB,REORD,AMPS,
     &           DYSEN,DUMMY,TITLE)

      CALL mma_deallocate(REORD)

      RETURN

100   FORMAT(1X,'There are contributions from several irreps, ',A,
     &       ' is not generated.')

      END
