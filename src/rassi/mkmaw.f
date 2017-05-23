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
      SUBROUTINE MKMAW(NLEV,NVERT,IDOWN,IDAW,IUP,IRAW,IMAW,
     &                 LTV,MIDLEV)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='MKMAW')
      DIMENSION IDOWN(NVERT,0:3),IDAW(NVERT,0:4)
      DIMENSION IUP(NVERT,0:3),IRAW(NVERT,0:4)
      DIMENSION IMAW(NVERT,0:3),LTV(-1:NLEV)
C Purpose: Construct the Modified Arc Weights



      CALL QENTER(ROUTINE)

C BEGIN BY CONSTRUCTING THE UPCHAIN TABLE IUP:
      DO IV=1,NVERT
        DO IC=0,3
          IUP(IV,IC)=0
        END DO
      END DO
      DO IU=1,NVERT-1
        DO IC=0,3
          IDWN=IDOWN(IU,IC)
          IF(IDWN.NE.0) IUP(IDWN,IC)=IU
        END DO
      END DO
C USE UPCHAIN TABLE TO CALCULATE THE REVERSE ARC WEIGHT TABLE:
      DO IC=0,3
        IRAW(1,IC)=0
      END DO
      IRAW(1,4)=1
      DO IV=2,NVERT
        ISUM=0
        DO IC=0,3
          IRAW(IV,IC)=0
          IU=IUP(IV,IC)
          IF(IU.NE.0) THEN
            IRAW(IV,IC)=ISUM
            ISUM=ISUM+IRAW(IU,4)
          END IF
        END DO
        IRAW(IV,4)=ISUM
      END DO
C USE IDAW,IRAW TABLES TO DETERMINE MIDLEV.
C ASSUMPTION IS THAT A BALANCED NUMBER OF UPPER/LOWER WALKS
C IS THE BEST CHOICE.
      MIDLEV=1
      MINW=1000000
      DO IL=1,NLEV-1
        NW=0
        DO IV=LTV(IL),LTV(IL-1)-1
          NW=NW+IRAW(IV,4)-IDAW(IV,4)
        END DO
        NW=ABS(NW)
        IF(NW.LT.MINW) THEN
          MIDLEV=IL
          MINW=NW
        END IF
      END DO
      MIDV1=LTV(MIDLEV)
      MIDV2=LTV(MIDLEV-1)-1
C COPY LOWER PART OF DIRECT ARC WEIGHT TABLE INTO IMAW:
      DO IV=MIDV1,NVERT
        DO IC=0,3
          IMAW(IV,IC)=IDAW(IV,IC)
        END DO
      END DO
C COPY UPPER PART OF REVERSE ARC WEIGHT TABLE INTO IMAW. HOWEVER,
C    NOTE THAT THE IMAW TABLE IS ACCESSED BY THE UPPER VERTEX.
      DO IU=1,MIDV1-1
        DO IC=0,3
          ID=IDOWN(IU,IC)
          IMAW(IU,IC)=0
          IF(ID.NE.0) IMAW(IU,IC)=IRAW(ID,IC)
        END DO
      END DO
C COPY UPPER PART OF REVERSE ARC WEIGHT TABLE INTO IMAW. HOWEVER,
C FINALLY, ADD AN OFFSET TO ARCS LEADING TO MIDLEVELS:
      ISUM=1
      DO IV=MIDV1,MIDV2
        DO IC=0,3
          IU=IUP(IV,IC)
          IF(IU.NE.0) IMAW(IU,IC)=ISUM+IMAW(IU,IC)
        END DO
        ISUM=ISUM+IRAW(IV,4)
      END DO
      DO IV=MIDV1,MIDV2
        DO IC=0,3
          IF(IDOWN(IV,IC).NE.0) IMAW(IV,IC)=ISUM+IMAW(IV,IC)
        END DO
        ISUM=ISUM+IDAW(IV,4)
      END DO

      CALL QEXIT(ROUTINE)
      RETURN
      END
