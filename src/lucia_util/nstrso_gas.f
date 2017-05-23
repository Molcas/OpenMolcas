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
* Copyright (C) 1994, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE NSTRSO_GAS(    NEL,  NORB1,  NORB2,  NORB3, NELMN1,
     &                       NELMX1, NELMN3, NELMX3,    IOC,   NORB,
     &                       NSTASO, ISTASO, NOCTYP,  NSMST,  IOTYP,
     &                        IPRNT)
*
* Number of strings per symmetry for group IOTYP
*
* Gas version, no check of type : set to 1
*
* Jeppe Olsen Winter of 1994
*
      IMPLICIT REAL*8           ( A-H,O-Z)
C     DIMENSION IOC(*),NSTASO(NOCTYP,NSMST)
      DIMENSION IOC(*),NSTASO(NSMST,*),ISTASO(NSMST,*)
*
      CALL ISETVC(NSTASO(1,IOTYP),0,NSMST)
      NTEST0 = 0
      NTEST = MAX(IPRNT,NTEST0)
      NSTRIN = 0
      IORB1F = 1
      IORB1L = IORB1F+NORB1-1
      IORB2F = IORB1L + 1
      IORB2L = IORB2F+NORB2-1
      IORB3F = IORB2L + 1
      IORB3L = IORB3F+NORB3-1
* Loop over possible partitionings between RAS1,RAS2,RAS3
      DO 1001 IEL1 = NELMX1,NELMN1,-1
      DO 1003 IEL3 = NELMN3,NELMX3, 1
       IF(IEL1.GT. NORB1 ) GOTO 1001
       IF(IEL3.GT. NORB3 ) GOTO 1003
       IEL2 = NEL - IEL1-IEL3
       IF(IEL2 .LT. 0 .OR. IEL2 .GT. NORB2 ) GOTO 1003
       IFRST1 = 1
* Loop over RAS 1 occupancies
  901  CONTINUE
         IF( IEL1 .NE. 0 ) THEN
           IF(IFRST1.EQ.1) THEN
            CALL ISTVC2(IOC(1),0,1,IEL1)
            IFRST1 = 0
           ELSE
             CALL NXTORD(IOC,IEL1,IORB1F,IORB1L,NONEW1)
             IF(NONEW1 .EQ. 1 ) GOTO 1003
           END IF
         END IF
         IF( NTEST .GE.500) THEN
           WRITE(6,*) ' RAS 1 string '
           CALL IWRTMA(IOC,1,IEL1,1,IEL1)
         END IF
         IFRST2 = 1
         IFRST3 = 1
* Loop over RAS 2 occupancies
  902    CONTINUE
           IF( IEL2 .NE. 0 ) THEN
             IF(IFRST2.EQ.1) THEN
              CALL ISTVC2(IOC(IEL1+1),IORB2F-1,1,IEL2)
              IFRST2 = 0
             ELSE
               CALL NXTORD(IOC(IEL1+1),IEL2,IORB2F,IORB2L,NONEW2)
               IF(NONEW2 .EQ. 1 ) THEN
                 IF(IEL1 .NE. 0 ) GOTO 901
                 IF(IEL1 .EQ. 0 ) GOTO 1003
               END IF
             END IF
           END IF
           IF( NTEST .GE.500) THEN
             WRITE(6,*) ' RAS 1 2 string '
             CALL IWRTMA(IOC,1,IEL1+IEL2,1,IEL1+IEL2)
           END IF
           IFRST3 = 1
* Loop over RAS 3 occupancies
  903      CONTINUE
             IF( IEL3 .NE. 0 ) THEN
               IF(IFRST3.EQ.1) THEN
                CALL ISTVC2(IOC(IEL1+IEL2+1),IORB3F-1,1,IEL3)
                IFRST3 = 0
               ELSE
                 CALL NXTORD(IOC(IEL1+IEL2+1),
     &           IEL3,IORB3F,IORB3L,NONEW3)
                 IF(NONEW3 .EQ. 1 ) THEN
                   IF(IEL2 .NE. 0 ) GOTO 902
                   IF(IEL1 .NE. 0 ) GOTO 901
                   GOTO 1003
                 END IF
               END IF
             END IF
             IF( NTEST .GE. 500) THEN
               WRITE(6,*) ' RAS 1 2 3 string '
               CALL IWRTMA(IOC,1,NEL,1,NEL)
             END IF
* Next string has been constructed , Enlist it !.
             NSTRIN = NSTRIN + 1
*. Symmetry of string
             ISYM = ISYMST(IOC,NEL)
C                   ISYMST(STRING,NEL)
*. occupation type of string
COLD         ITYP = IOCTP2(IOC,NEL,IOTYP)
C                   IOCTP2(STRING,NEL)
*
             NSTASO(ISYM,IOTYP) = NSTASO(ISYM,IOTYP)+ 1
*
           IF( IEL3 .NE. 0 ) GOTO 903
           IF( IEL3 .EQ. 0 .AND. IEL2 .NE. 0 ) GOTO 902
           IF( IEL3 .EQ. 0 .AND. IEL2 .EQ. 0 .AND. IEL1 .NE. 0)
     &     GOTO 901
 1003 CONTINUE
 1001 CONTINUE
*
*. The corresponding offset
*
      DO ISM = 1, NSMST
        IF(ISM.EQ.1) THEN
          ISTASO(ISM,IOTYP) = 1
        ELSE
          ISTASO(ISM,IOTYP) = ISTASO(ISM-1,IOTYP)+NSTASO(ISM-1,IOTYP)
        END IF
      END DO

      IF(NTEST.GE.5)
     &WRITE(6,*) ' Number of strings generated   ', NSTRIN
      IF(NTEST .GE. 10 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' Number of strings per sym for group = ', IOTYP
        WRITE(6,*) '================================================'
        CALL IWRTMA(NSTASO(1,IOTYP),1,NSMST,1,NSMST)
        WRITE(6,*) ' Offset for given symmetry for group = ', IOTYP
        WRITE(6,*) '================================================'
        CALL IWRTMA(ISTASO(1,IOTYP),1,NSMST,1,NSMST)
      END IF
C
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer(NORB)
        CALL Unused_integer(NOCTYP)
      END IF
      END
