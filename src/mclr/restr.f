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
      SUBROUTINE RESTR_MCLR(NVERT0,IDRT0,IDOWN0,IVER,
     &                 LV1RAS,LV3RAS,LM1RAS,LM3RAS,NVERT)
C
C     PURPOSE: PUT THE RAS CONSTRAINT TO THE DRT TABLE BY
C              CREATING A MASK
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION IDRT0(NVERT0,5),IDOWN0(NVERT0,0:3),IVER(NVERT0)
C
      PARAMETER (LTAB=1,NTAB=2)
      DIMENSION IIOR(0:3,0:3)
      DIMENSION IIAND(0:3,0:3)
      DATA IIOR  / 0,1,2,3,1,1,3,3,2,3,2,3,3,3,3,3 /
      DATA IIAND / 0,0,0,0,0,1,0,1,0,0,2,2,0,1,2,3 /
C
C
C     LOOP OVER ALL VERTICES AND CHECK ON RAS CONDITIONS
C     CREATE MASK
C
      DO 10 IV=1,NVERT0
        LEV=IDRT0(IV,LTAB)
        N=IDRT0(IV,NTAB)
        IVER(IV)=0
        IF((LEV.EQ.LV1RAS).AND.(N.GE.LM1RAS)) IVER(IV)=1
        IF((LEV.EQ.LV3RAS).AND.(N.GE.LM3RAS)) IVER(IV)=IVER(IV)+2
10    CONTINUE
C
C     NOW LOOP FORWARDS, MARKING THOSE VERTICES CONNECTED FROM ABOVE.
C     SINCE IVER WAS INITIALIZED TO ZERO, NO CHECKING IS NEEDED.
C
      DO 20 IV=1,NVERT0-1
        IVV=IVER(IV)
        DO 20 IC=0,3
          ID=IDOWN0(IV,IC)
          IF(ID.EQ.0) GOTO 20
          IVER(ID)=IIOR(IVER(ID),IVV)
20    CONTINUE
C
C     THEN LOOP BACKWARDS. SAME RULES, EXCEPT THAT CONNECTIVITY
C     SHOULD BE PROPAGATED ONLY ABOVE THE RESTRICTION LEVELS.
C
      DO 40 IV=NVERT0-1,1,-1
        LEV=IDRT0(IV,LTAB)
        MASK=0
        IF(LEV.GT.LV1RAS) MASK=1
        IF(LEV.GT.LV3RAS) MASK=MASK+2
        IVV=IVER(IV)
        DO 30 IC=0,3
          ID=IDOWN0(IV,IC)
          IF(ID.EQ.0) GOTO 30
          IVD=IVER(ID)
          IVV=IIOR(IVV,IIAND(MASK,IVD))
30      CONTINUE
        IVER(IV)=IVV
40    CONTINUE
C
C     WE ARE NOW INTERESTED ONLY IN VERTICES CONNECTED BOTH TO
C     ALLOWED VERTICES FOR RAS-SPACE 1 AND RAS-SPACE 3.
C     THOSE ARE NUMBERED IN ASCENDING ORDER, THE REST ARE ZEROED.
C
      NVERT=0
      DO 50 IV=1,NVERT0
        IF(IVER(IV).EQ.3) THEN
          NVERT=NVERT+1
          IVER(IV)=NVERT
        ELSE
          IVER(IV)=0
        ENDIF
50    CONTINUE
C
C
C     EXIT
C
      RETURN
      END
