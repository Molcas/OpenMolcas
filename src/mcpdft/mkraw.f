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
      SUBROUTINE MKRAW_m(IDOWN,IUP,IRAW,IPRINT)
C
C     PURPOSE: CONSTRUCT UPCHAIN INDEX TABLE AND REVERSE ARC WEIGHTS
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "gugx.fh"
#include "rasdim.fh"
#include "general.fh"
#include "output_ras.fh"
      Parameter(Routine='MKRAW')
C
      DIMENSION IDOWN(NVERT,0:3),IUP(NVERT,0:3),IRAW(NVERT,0:4)

C
C     BEGIN BY CONSTRUCTING THE UPCHAIN TABLE IUP:
C
      DO IV=1,NVERT
        DO IC=0,3
          IUP(IV,IC)=0
        END DO
      END DO
      DO IU=1,NVERT-1
        DO IC=0,3
          IDWN=IDOWN(IU,IC)
          IF(IDWN.EQ.0) GOTO 20
          IUP(IDWN,IC)=IU
20        CONTINUE
        END DO
      END DO
C
      IF( IPRINT.GE.5 ) THEN
        Write(LF,*)
        Write(LF,*)' THE UPCHAIN TABLE IN MKRAW:'
        DO IV=1,NVERT
          Write(LF,'(1X,I4,5X,4(1X,I6))') IV,(IUP(IV,IC),IC=0,3)
        END DO
        Write(LF,*)
      ENDIF
C
C     USE UPCHAIN TABLE TO CALCULATE THE REVERSE ARC WEIGHT TABLE:
C
      DO IC=0,3
        IRAW(1,IC)=0
      END DO
      IRAW(1,4)=1
      DO IV=2,NVERT
        ISUM=0
        DO IC=0,3
          IRAW(IV,IC)=0
          IU=IUP(IV,IC)
          IF(IU.EQ.0) GOTO 120
          IRAW(IV,IC)=ISUM
          ISUM=ISUM+IRAW(IU,4)
120     CONTINUE
        END DO
        IRAW(IV,4)=ISUM
      END DO
C
      IF( IPRINT.GE.5 ) THEN
        Write(LF,*)
        Write(LF,*)' THE REVERSE ARC WEIGHT TABLE IN MKRAW:'
        DO IV=1,NVERT
          Write(LF,'(1X,I4,5X,5(1X,I6))') IV,(IRAW(IV,IC),IC=0,4)
        END DO
        Write(LF,*)
      ENDIF

      RETURN
      END
