!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!#define _DEBUGPRINT_
      SUBROUTINE MKMID(NVERT,NLEV,IDRT,IDAW,IRAW,LTV,      &
                       MIDLEV, NMIDV, MIDV1, MIDV2, MXUP, MXDWN)
!     PURPOSE: FIND THE MIDLEVEL
!
#ifdef _DEBUGPRINT_
      use Definitions, only: LF => u6
#endif
      IMPLICIT None
!
      Integer NVERT, NLEV, MIDLEV, NMIDV, MIDV1, MIDV2, MXUP, MXDWN
      Integer IDRT(NVERT,5)
      Integer IDAW(NVERT,0:4)
      Integer IRAW(NVERT,0:4)
      Integer LTV(-1:NLEV)

      Integer, Parameter:: LTAB=1
      Integer IV, LEV, MINW, MV, NW, IL
!
!     SET UP A LEVEL-TO-VERTEX TABLE, LTV, AND IDENTIFY MIDVERTICES:
!
      LTV(:)=0
!
      DO IV=1,NVERT
        LEV=IDRT(IV,LTAB)
        LTV(LEV)=LTV(LEV)+1
      End Do
!
      DO LEV=NLEV,0,-1
        LTV(LEV-1)=LTV(LEV-1)+LTV(LEV)
      End Do
!
      DO LEV=-1,NLEV-1
        LTV(LEV)=1+LTV(LEV+1)
      End Do
!
!     USE IDAW,IRAW TABLES TO DETERMINE MIDLEV.
!     THE ASSUMPTION IS THAT A BALANCED NUMBER OF UPPER/LOWER WALKS
!     IS THE BEST CHOICE.
!
!hrl 980529 fix for nLev=0 (no orbitals in any active space)
!     Since LTV(-1:nLev)  and the statement after the loop
!     MidV1=LTV(MidLev) we have the condition MidLev>=nLev
!     Hence MidLev=1 is inappropriate for nLev=0
!     MIDLEV=1
!
      If (nLev==0) Then
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
        IF(NW.GE.MINW) Cycle
        MIDLEV=IL
        MINW=NW
      END DO
      MIDV1=LTV(MIDLEV)
      MIDV2=LTV(MIDLEV-1)-1
      NMIDV=1+MIDV2-MIDV1
!
!     NOW FIND THE MAX NUMBERS OF UPPER AND LOWER WALKS. RESPECTIVELY
!     (DISREGARDING SYMMETRY)
!
      MXUP=0
      MXDWN=0
      DO MV=MIDV1,MIDV2
        if(MXUP<IRAW(MV,4)) MXUP=IRAW(MV,4)
        if(MXDWN<IDAW(MV,4)) MXDWN=IDAW(MV,4)
      END DO
!
#ifdef _DEBUGPRINT_
      Write(LF,*)
      Write(LF,'(A,I3)')' MIDLEVEL =             ',MIDLEV
      Write(LF,'(A,I3)')' NUMBER OF MIDVERTICES =',NMIDV
      Write(LF,'(A,I3)')' FIRST MIDVERTEX =      ',MIDV1
      Write(LF,'(A,I3)')' LAST MIDVERTEX =       ',MIDV2
      Write(LF,'(A,I3)')' MAX. NO UPPER WALKS=   ',MXUP
      Write(LF,'(A,I3)')' MAX. NO LOWER WALKS=   ',MXDWN
      Write(LF,*)
#endif

      END SUBROUTINE MKMID
