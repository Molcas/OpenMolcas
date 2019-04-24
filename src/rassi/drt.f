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
      SUBROUTINE DRT_RASSI(NVERT0,IDRT0,IDOWN0,NWVER,NVERT,IDRT,IDOWN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IDRT0(NVERT0,5),IDOWN0(NVERT0,0:3)
      DIMENSION NWVER(NVERT0)
      DIMENSION IDRT(NVERT,5),IDOWN(NVERT,0:3)
C PURPOSE: IDRT0 AND IDOWN0 ARE PALDUS DRT TABLES FOR THE FULL
C    GRAPH. NWVER IS A REINDICING TABLE, WHERE NWVER(IV)=0 IF VERTEX
C    IV SHOULD BE EXCLUDED, AND ELSE IS THE NEW INDEX FOR THAT VERTEX.
C    OUTPUT IS IDRT AND IDOWN, THE TABLES FOR THE RESTRICTED GRAPH.
      DO IV=1,NVERT0
        IVNEW=NWVER(IV)
        IF(IVNEW.EQ.0) GOTO 30
        DO ITAB=1,5
          IDRT(IVNEW,ITAB)=IDRT0(IV,ITAB)
        END DO
        DO IC=0,3
          ID=IDOWN0(IV,IC)
          IDNEW=0
          IF(ID.NE.0) IDNEW=NWVER(ID)
          IDOWN(IVNEW,IC)=IDNEW
        END DO
30    CONTINUE
      END DO
      RETURN
      END
