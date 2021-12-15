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
      SUBROUTINE RMVERT(NLEV,NVERT,IDRT,IDOWN,NLIM,NWVERT)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='RMVERT')
#include "WrkSpc.fh"
      DIMENSION IDRT(NVERT,5),IDOWN(NVERT,0:3)
      DIMENSION NWVERT(NVERT)
      PARAMETER (LTAB=1,NTAB=2)
      DIMENSION NLIM(NLEV)
C Purpose: Remove vertices from a DRT table.




      CALL GETMEM('Conn','Allo','Inte',LCONN,NVERT)
C KILL VERTICES THAT DO NOT OBEY RESTRICTIONS.
      DO IV=1,NVERT-1
        NWVERT(IV)=1
        L=IDRT(IV,LTAB)
        N=IDRT(IV,NTAB)
        IF(N.LT.NLIM(L)) NWVERT(IV)=0
      END DO
      NWVERT(NVERT)=1

C     OPEN(UNIT=10,STATUS='NEW',FILE='SGraph.ps')
C     call DRAWSG(10,NLEV,NVERT,IDRT,IDOWN,NWVERT)

  10  CONTINUE
C REMOVE ARCS HAVING A DEAD UPPER OR LOWER VERTEX.
C COUNT THE NUMBER OF ARCS REMOVED OR VERTICES KILLED.
      NCHANGES=0
      DO IV=1,NVERT-1
        IF(NWVERT(IV).EQ.0) THEN
          DO IC=0,3
            ID=IDOWN(IV,IC)
            IF(ID.GT.0) THEN
              IDOWN(IV,IC)=0
              NCHANGES=NCHANGES+1
            END IF
          END DO
        ELSE
          NLD=0
          DO IC=0,3
            ID=IDOWN(IV,IC)
            IF(ID.GT.0) THEN
              IF(NWVERT(ID).EQ.0) THEN
                IDOWN(IV,IC)=0
                NCHANGES=NCHANGES+1
              ELSE
                NLD=NLD+1
              END IF
            END IF
          END DO
          IF(NLD.EQ.0) THEN
            NWVERT(IV)=0
            NCHANGES=NCHANGES+1
          END IF
        END IF
      END DO
C ALSO CHECK ON CONNECTIONS FROM ABOVE:
      IWORK(LCONN)=NWVERT(1)
      DO IV=2,NVERT
        IWORK(LCONN-1+IV)=0
      END DO
      DO IV=1,NVERT-1
        IF(NWVERT(IV).EQ.1) THEN
          DO IC=0,3
            ID=IDOWN(IV,IC)
            IF(ID.GT.0 .AND. NWVERT(ID).EQ.1) THEN
                IWORK(LCONN-1+ID)=1
            END IF
          END DO
        END IF
      END DO
      DO IV=1,NVERT
        IF(NWVERT(IV).EQ.1 .AND. IWORK(LCONN-1+IV).EQ.0) THEN
          NWVERT(IV)=0
          NCHANGES=NCHANGES+1
        END IF
      END DO
C     call DRAWSG(10,NLEV,NVERT,IDRT,IDOWN,NWVERT)
C ANY CHANGES? RERUN.
      IF(NCHANGES.GT.0) GOTO 10
      CALL GETMEM('Conn','Free','Inte',LCONN,NVERT)
C     Close(10)

C IF NO CHANGES, THE REMAINING GRAPH IS VALID.
C EVERY VERTEX OBEYS THE RESTRICTIONS. EVERY VERTEX IS
C CONNECTED ABOVE AND BELOW (EXCEPTING THE TOP AND BOTTOM)
C TO OTHER CONFORMING VERTICES.
C THE PROCEDURE IS GUARANTEED TO FIND A STABLE SOLUTIONS,
C SINCE EACH ITERATION REMOVES ARCS AND/OR VERTICES FROM THE
C FINITE NUMBER WE STARTED WITH.
      IF(NWVERT(1).EQ.0) THEN
        WRITE(6,*)'RASSI/RMVERT: Too severe restrictions.'
        WRITE(6,*)'Not one single configuration is left.'
        CALL ABEND()
      END IF

      NV=0
      DO IV=1,NVERT
        IF(NWVERT(IV).EQ.1) THEN
          NV=NV+1
          NWVERT(IV)=NV
        END IF
      END DO
      NVERT=NV


      RETURN
      END
