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
* Copyright (C) 1995,2000, Jeppe Olsen                                 *
************************************************************************
      SUBROUTINE ORBINH1(IORBINH1,IORBINH1_NOCCSYM,NTOOBS,NTOOB,NSMOB)
*
* Obtain array of 2 orbital indeces,
* for symmetry packed matrices
*
* IORBINH1 : Lower half packed
* IORBINH1_NOCCSYM : Complete blocks
*
* resulting indeces are with respect to start of given symmetry block
* while input orbital indeces are absolute and in symmetry order
*
* Jeppe Olsen, March 1995
*              ORBINH1_NOCCSYM added August 2000
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION NTOOBS(NSMOB)
*. output
      DIMENSION IORBINH1(NTOOB,NTOOB), IORBINH1_NOCCSYM(NTOOB,NTOOB)
*
C?    WRITE(6,*) ' ORBINH1 speaking '
C?    WRITE(6,*) ' NSMOB NTOOB ',NSMOB,NTOOB
C?    WRITE(6,*) ' NTOOBS '
C?    CALL IWRTMA(NTOOBS,1,NSMOB,1,NSMOB)
*. To eliminate annoying and incorrect compiler warnings
      IOFF = 0
      JOFF = 0
      INDEX = 0

*. Loop over symmetries of orbitals

      DO ISM = 1, NSMOB
        IF(ISM.EQ.1) THEN
          IOFF = 1
        ELSE
          IOFF = IOFF + NTOOBS(ISM-1)
        END IF
        DO JSM = 1, NSMOB
          IF(JSM.EQ.1) THEN
            JOFF = 1
          ELSE
            JOFF = JOFF + NTOOBS(JSM-1)
          END IF
C?        WRITE(6,*) ' ISM JSM IOFF JOFF', ISM,JSM,IOFF,JOFF
          DO IORB = 1, NTOOBS(ISM)
            IABS = IOFF -1 + IORB
            DO JORB = 1, NTOOBS(JSM)
              JABS = JOFF -1 + JORB
C?            write(6,*) ' IORB JORB IABS JABS ',IORB,JORB,IABS,JABS
              IF(ISM.GT.JSM) THEN
                INDEX = (IORB-1)*NTOOBS(JSM) + JORB
              ELSE IF(ISM.EQ.JSM) THEN
                IF(IORB.GE.JORB) THEN
                  INDEX = IORB*(IORB-1)/2 + JORB
                ELSE
                  INDEX = JORB*(JORB-1)/2 + IORB
                END IF
              ELSE IF(ISM.LT.JSM) THEN
                INDEX = (JORB-1)*NTOOBS(ISM) + IORB
              END IF
              INDEX_NOCCSYM = (IORB-1)*NTOOBS(JSM) + JORB
              IORBINH1(IABS,JABS) = INDEX
              IORBINH1_NOCCSYM(IABS,JABS) = INDEX_NOCCSYM
            END DO
          END DO
*. End of loops over orbital indeces
        END DO
      END DO
*. End of loop over orbital symmetries
*
      NTEST = 000
      IF(NTEST .GE. 100 ) THEN
        WRITE(6,*) ' IORBINH1 matrix delivered from ORBINH1'
        CALL IWRTMA(IORBINH1,NTOOB,NTOOB,NTOOB,NTOOB)
      END IF
*
      RETURN
      END
