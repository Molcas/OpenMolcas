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
      SUBROUTINE MKMID(SGS)
!     PURPOSE: FIND THE MIDLEVEL
!
#ifdef _DEBUGPRINT_
      use Definitions, only: LF => u6
#endif
      use struct, only: SGStruct
      IMPLICIT None
!
      Type (SGStruct) SGS

      Integer IV, MINW, MV, NW, IL
#ifdef _DEBUGPRINT_
      Integer NMIDV
#endif

      Associate (nVert=>SGS%nVert, nLev=>SGS%nLev, MidLev => SGS%MidLev, &
                 MvSta=>SGS%MVSta, MVEnd=>SGS%MVEnd,   &
                 MxUp=>SGS%MxUP, MxDWN=>SGS%MxDwn, iDAW=>SGS%Daw,        &
                 iRaw=>SGS%Raw, LTV=>SGS%LTV)
!
!     USE IDAW,IRAW TABLES TO DETERMINE MIDLEV.
!     THE ASSUMPTION IS THAT A BALANCED NUMBER OF UPPER/LOWER WALKS
!     IS THE BEST CHOICE.
!
!hrl 980529 fix for nLev=0 (no orbitals in any active space)
!     Since LTV(-1:nLev)  and the statement after the loop
!     MVSta=LTV(MidLev) we have the condition MidLev>=nLev
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
      MVSta=LTV(MIDLEV)
      MVEnd=LTV(MIDLEV-1)-1
!
!     NOW FIND THE MAX NUMBERS OF UPPER AND LOWER WALKS. RESPECTIVELY
!     (DISREGARDING SYMMETRY)
!
      MXUP=0
      MXDWN=0
      DO MV=MVSta,MVEnd
        if(MXUP<IRAW(MV,4)) MXUP=IRAW(MV,4)
        if(MXDWN<IDAW(MV,4)) MXDWN=IDAW(MV,4)
      END DO
!
#ifdef _DEBUGPRINT_
      NMIDV=>MVEnd-MVSta+1
      Write(LF,*)
      Write(LF,'(A,I3)')' MIDLEVEL =             ',MIDLEV
      Write(LF,'(A,I3)')' NUMBER OF MIDVERTICES =',NMIDV
      Write(LF,'(A,I3)')' FIRST MIDVERTEX =      ',MVSta
      Write(LF,'(A,I3)')' LAST MIDVERTEX =       ',MVEnd
      Write(LF,'(A,I3)')' MAX. NO UPPER WALKS=   ',MXUP
      Write(LF,'(A,I3)')' MAX. NO LOWER WALKS=   ',MXDWN
      Write(LF,*)
      d Associate
#endif

      End Associate

      END SUBROUTINE MKMID
