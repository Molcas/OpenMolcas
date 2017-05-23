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
      SUBROUTINE MKMID(NVERT,NLEV,IDRT,IDOWN,IDAW,IUP,IRAW,LTV,
     &                 MIDLEV,NMIDV,MIDV1,MIDV2,MXUP,MXDWN,IPRINT)
C
C     PURPOSE: FIND THE MIDLEVEL
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION IDRT(NVERT,5)
      DIMENSION IDOWN(NVERT,0:3),IDAW(NVERT,0:4)
      DIMENSION IUP(NVERT,0:3),IRAW(NVERT,0:4)
      DIMENSION LTV(-1:NLEV)
C
C
C     SET UP A LEVEL-TO-VERTEX TABLE, LTV, AND IDENTIFY MIDVERTICES:
C
      LTAB=1
      DO  LEV=-1,NLEV
        LTV(LEV)=0
      END DO
      DO IV=1,NVERT
        LEV=IDRT(IV,LTAB)
        LTV(LEV)=LTV(LEV)+1
      END DO
      DO  LEV=NLEV,0,-1
       LTV(LEV-1)=LTV(LEV-1)+LTV(LEV)
      END DO
      DO LEV=-1,NLEV-1
        LTV(LEV)=1+LTV(LEV+1)
      END DO
C
C     USE IDAW,IRAW TABLES TO DETERMINE MIDLEV.
C     THE ASSUMPTION IS THAT A BALANCED NUMBER OF UPPER/LOWER WALKS
C     IS THE BEST CHOICE.
C
      IF (NLEV.eq.0) Then
       MIDLEV=0
      Else
      MIDLEV=1
      End If
      MINW=1000000
      DO 160 IL=1,NLEV-1
        NW=0
        DO 150 IV=LTV(IL),LTV(IL-1)-1
          NW=NW+IRAW(IV,4)-IDAW(IV,4)
150     CONTINUE
        NW=ABS(NW)
        IF(NW.GE.MINW) GOTO 160
        MIDLEV=IL
        MINW=NW
160   CONTINUE
      MIDV1=LTV(MIDLEV)
      MIDV2=LTV(MIDLEV-1)-1
      NMIDV=1+MIDV2-MIDV1
C
C     NOW FIND THE MAX NUMBERS OF UPPER AND LOWER WALKS. RESPECTIVELY
C     (DISREGARDING SYMMETRY)
C
      MXUP=0
      MXDWN=0
      DO 200 MV=MIDV1,MIDV2
        MXUP=MAX(MXUP,IRAW(MV,4))
        MXDWN=MAX(MXDWN,IDAW(MV,4))
200   CONTINUE
C
      IF( IPRINT.GE.5 ) THEN
        WRITE(6,*)
        WRITE(6,'(A,I3)')' MIDLEVEL =             ',MIDLEV
        WRITE(6,'(A,I3)')' NUMBER OF MIDVERTICES =',NMIDV
        WRITE(6,'(A,I3)')' FIRST MIDVERTEX =      ',MIDV1
        WRITE(6,'(A,I3)')' LAST MIDVERTEX =       ',MIDV2
        WRITE(6,'(A,I3)')' MAX. NO UPPER WALKS=   ',MXUP
        WRITE(6,'(A,I3)')' MAX. NO LOWER WALKS=   ',MXDWN
        WRITE(6,*)
      ENDIF
C
C
C     EXIT
C
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer_array(IDOWN)
        CALL Unused_integer_array(IUP)
      END IF
      END
