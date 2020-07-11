************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE MKRAW_MCLR(NVERT,IDOWN,IDAW,IUP,IRAW,IPRINT)
C
C     PURPOSE: CONSTRUCT UPCHAIN INDEX TABLE AND REVERSE ARC WEIGHTS
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION IDOWN(NVERT,0:3),IDAW(NVERT,0:4)
      DIMENSION IUP(NVERT,0:3),IRAW(NVERT,0:4)
C
C
C     BEGIN BY CONSTRUCTING THE UPCHAIN TABLE IUP:
C
      DO 10 IV=1,NVERT
        DO 11 IC=0,3
          IUP(IV,IC)=0
11      CONTINUE
10    CONTINUE
      DO 30 IU=1,NVERT-1
        DO 20 IC=0,3
          IDWN=IDOWN(IU,IC)
          IF(IDWN.EQ.0) GOTO 20
          IUP(IDWN,IC)=IU
20      CONTINUE
30    CONTINUE
C
      IF( IPRINT.GE.5 ) THEN
        WRITE(6,*)
        WRITE(6,*)' THE UPCHAIN TABLE IN MKRAW_MCLR:'
        DO 40 IV=1,NVERT
          WRITE(6,'(1X,I4,5X,4(1X,I6))') IV,(IUP(IV,IC),IC=0,3)
40      CONTINUE
        WRITE(6,*)
      ENDIF
C
C     USE UPCHAIN TABLE TO CALCULATE THE REVERSE ARC WEIGHT TABLE:
C
      DO 110 IC=0,3
        IRAW(1,IC)=0
110   CONTINUE
      IRAW(1,4)=1
      DO 130 IV=2,NVERT
        ISUM=0
        DO 120 IC=0,3
          IRAW(IV,IC)=0
          IU=IUP(IV,IC)
          IF(IU.EQ.0) GOTO 120
          IRAW(IV,IC)=ISUM
          ISUM=ISUM+IRAW(IU,4)
120     CONTINUE
        IRAW(IV,4)=ISUM
130   CONTINUE
C
      IF( IPRINT.GE.5 ) THEN
        WRITE(6,*)
        WRITE(6,*)' THE REVERSE ARC WEIGHT TABLE IN MKRAW_MCLR:'
        DO 140 IV=1,NVERT
          WRITE(6,'(1X,I4,5X,5(1X,I6))') IV,(IRAW(IV,IC),IC=0,4)
140     CONTINUE
        WRITE(6,*)
      ENDIF
C
C
C     EXIT
C
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer_array(IDAW)
      END
