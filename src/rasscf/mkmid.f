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
      SUBROUTINE MKMID(IDRT,IDAW,IRAW,LTV,IPRINT)
C     PURPOSE: FIND THE MIDLEVEL
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "gugx.fh"
#include "rasdim.fh"
#include "general.fh"
#include "output_ras.fh"
      Parameter(Routine='MKMID')
C
      DIMENSION IDRT(NVERT,5)
      DIMENSION IDAW(NVERT,0:4)
      DIMENSION IRAW(NVERT,0:4)
      DIMENSION LTV(-1:NLEV)

      Call qEnter(Routine)

C
C     SET UP A LEVEL-TO-VERTEX TABLE, LTV, AND IDENTIFY MIDVERTICES:
C
      LTAB=1
      DO LEV=-1,NLEV
        LTV(LEV)=0
      End Do
*
      DO IV=1,NVERT
        LEV=IDRT(IV,LTAB)
        LTV(LEV)=LTV(LEV)+1
      End Do
*
      DO LEV=NLEV,0,-1
        LTV(LEV-1)=LTV(LEV-1)+LTV(LEV)
      End Do
*
      DO LEV=-1,NLEV-1
        LTV(LEV)=1+LTV(LEV+1)
      End Do
C
C     USE IDAW,IRAW TABLES TO DETERMINE MIDLEV.
C     THE ASSUMPTION IS THAT A BALANCED NUMBER OF UPPER/LOWER WALKS
C     IS THE BEST CHOICE.
C
*hrl 980529 fix for nLev=0 (no orbitals in any active space)
*     Since LTV(-1:nLev)  and the statement after the loop
*     MidV1=LTV(MidLev) we have the condition MidLev>=nLev
*     Hence MidLev=1 is inappropriate for nLev=0
*     MIDLEV=1
*
      If (nLev.eq.0) Then
         MIDLEV=0
      Else
         MIDLEV=1
      End If
      MINW=1000000
      DO IL=1,NLEV-1
        NW=0
        DO IV=LTV(IL),LTV(IL-1)-1
          NW=NW+IRAW(IV,4)-IDAW(IV,4)
        END DO
        NW=ABS(NW)
        IF(NW.GE.MINW) GOTO 160
        MIDLEV=IL
        MINW=NW
160   CONTINUE
      END DO
      MIDV1=LTV(MIDLEV)
      MIDV2=LTV(MIDLEV-1)-1
      NMIDV=1+MIDV2-MIDV1
C
C     NOW FIND THE MAX NUMBERS OF UPPER AND LOWER WALKS. RESPECTIVELY
C     (DISREGARDING SYMMETRY)
C
      MXUP=0
      MXDWN=0
      DO MV=MIDV1,MIDV2
c        MXUP=MAX(MXUP,IRAW(MV,4))
        if(MXUP.lt.IRAW(MV,4)) MXUP=IRAW(MV,4)
c        MXDWN=MAX(MXDWN,IDAW(MV,4))
        if(MXDWN.lt.IDAW(MV,4)) MXDWN=IDAW(MV,4)
      END DO
C
      IF( IPRINT.GE.5 ) THEN
        Write(LF,*)
        Write(LF,'(A,I3)')' MIDLEVEL =             ',MIDLEV
        Write(LF,'(A,I3)')' NUMBER OF MIDVERTICES =',NMIDV
        Write(LF,'(A,I3)')' FIRST MIDVERTEX =      ',MIDV1
        Write(LF,'(A,I3)')' LAST MIDVERTEX =       ',MIDV2
        Write(LF,'(A,I3)')' MAX. NO UPPER WALKS=   ',MXUP
        Write(LF,'(A,I3)')' MAX. NO LOWER WALKS=   ',MXDWN
        Write(LF,*)
      ENDIF

      Call qExit(Routine)
      RETURN
      END
