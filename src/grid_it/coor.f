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
      Subroutine MyCoor(isAuto,
     &      Ox,Oy,Oz, Rx,Ry,Rz, iGx,iGy,iGz, iMagic,iCustOrig)
************************************************************************
* Adapted from SAGIT to work with OpenMolcas (October 2020)            *
************************************************************************
*                                                                      *
*   Read Coordinates and calculate a cub for grid                      *
*   isAuto=1 - real job, else only print                               *
*   Origin(3) fix cub in space                                         *
*   Rx,Ry,Rz - size of cub                                             *
*   iGx,iGy,iGz - net                                                  *
*   iMagic = magic guess for net                                       *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "WrkSpc.fh"
#include "grid.fh"
#include "real.fh"
c      Dimension iOper(8)
c      Dimension RotVec(3)
      Character*(LENIN) Byte4
c      Character*4 tt
      Character*128 Line
cVV: the constant is used in all GV packages
#define R529 0.52917721067d0
*----------------------------------------------------------------------*
*     Prologue                                                         *
*----------------------------------------------------------------------*
       if (isLUSCUS.eq.1) then
         x529=R529
       else
         x529=One
       endif
       lw2=1
      Call get_iScalar('nSym',nSym)
      Call Get_nAtoms_All(nAtoms)
      Call Get_LblCnt_All(AtomLbl)
      Call GetMem('Coor','ALLO','REAL',ipCoor,3*nAtoms)
      Call Get_Coord_All(Work(ipCoor),nAtoms)
      nCenter=nAtoms
* Check : are any strange names?
      NoOrig=0
      Do iAt=0,nCenter-1
       write(line,'(a)') AtomLbl(lw2+iAt)
       if(index( line,'Ori').ne.0) NoOrig=NoOrig+1
      enddo
*
      IF (ISLUSCUS .EQ. 1 .AND. ISBINARY .EQ. 3) THEN
        WRITE(LINE,'(2X,I8)') NCENTER-NOORIG
        CALL PRINTLINE(LID, LINE, 10,0)
        CALL PRINTLINE(LID, LINE, 0,0)
        if(isUHF.eq.1) then
        CALL PRINTLINE(LID_ab, LINE, 10,0)
        CALL PRINTLINE(LID_ab, LINE, 0,0)
        endif
        DO IAT=0,NCENTER-1
          WRITE(LINE,'(A)') ATOMLBL(LW2+IAT)
          IF(INDEX( LINE,'ORI').EQ.0) THEN
          Byte4=ATOMLBL(LW2+IAT)(1:2)
          if(index ('0123456789',Byte4(2:2)).ne.0) Byte4(2:2)=' '
            WRITE(LINE,'(1X,A2,2X,3F15.8)') Byte4,
     &      WORK(IPCOOR+3*IAT)*x529,WORK(IPCOOR+3*IAT+1)*x529,
     &      WORK(IPCOOR+3*IAT+2)*x529
            CALL PRINTLINE(LID, LINE, 50,0)
            if(isUHF.eq.1) then
            CALL PRINTLINE(LID_ab, LINE, 50,0)
            endif
          ENDIF
        ENDDO
      END IF

      Write(Line,'(A,I8)') 'Natom= ',nCenter-NoOrig
      If (isBinary.eq.1) Then
        Write(LuVal) Trim(Line)
        If (isUHF.eq.1) Write(LuVal_ab) Line(1:15)
      Else
        Write(LuVal,'(A)') Trim(Line)
        If (isUHF.eq.1) Write(LuVal_ab,'(A)') Line(1:15)
      End If
      Do iAt=0,nCenter-1
        If (Index(Line,'Ori').eq.0) Then
          Write(line,'(A,2X,3F15.8)') AtomLbl(lw2+iAt),
     &       Work(ipCoor+3*iAt  ),
     &       Work(ipCoor+3*iAt+1),
     &       Work(ipCoor+3*iAt+2)
        Else
          Write(Line,'(A)') AtomLbl(lw2+iAt)
        End If
        If (isBinary.eq.1) Then
          Write(LuVal) Trim(Line)
          If(isUHF.eq.1) Write(LuVal_ab) Trim(Line)
        Else
          Write(LuVal,'(A)') Trim(Line)
          If(isUHF.eq.1) Write(LuVal_ab,'(A)') Trim(Line)
        End If
      End Do

      if(iCustOrig.eq.1) goto 500

      if(isAuto.eq.1) Then
*----------------------------------------------------------------------*
*     Find Cub parameters                                              *
*                 Ox->RxMin, Rx->RxMax
*----------------------------------------------------------------------*
      Ox=9999.9D0
      Oy=9999.9D0
      Oz=9999.9D0
      Rx=-9999.9D0
      Ry=-9999.9D0
      Rz=-9999.9D0

      Do iAt=0,nCenter-1
       rrx=Work(ipCoor+3*iAt)
       rry=Work(ipCoor+3*iAt+1)
       rrz=Work(ipCoor+3*iAt+2)
       if(rrx.lt.Ox) Ox=rrx
       if(rrx.gt.Rx) Rx=rrx
       if(rry.lt.Oy) Oy=rry
       if(rry.gt.Ry) Ry=rry
       if(rrz.lt.Oz) Oz=rrz
       if(rrz.gt.Rz) Rz=rrz
      End Do
       Rx=Rx-Ox
       Ry=Ry-Oy
       Rz=Rz-Oz
*
* and now, expand this cub to place all atoms inside
*
      Ox=dble(int(Ox-TheGap))
      Oy=dble(int(Oy-TheGap))
      Oz=dble(int(Oz-TheGap))
      Rx=dble(int(Rx+Two*TheGap))
      Ry=dble(int(Ry+Two*TheGap))
      Rz=dble(int(Rz+Two*TheGap))
* finish of iAuto.
      endif
*----------------------------------------------------------------------*
*     Calculate corrected coords
*----------------------------------------------------------------------*
*
* make a stupid Patch: Cerius2 works well only with even nets!
*
       Rx=dble(int(Rx)/2*2)
       Ry=dble(int(Ry)/2*2)
       Rz=dble(int(Rz)/2*2)

       if(iMagic.gt.0) then
           iGx=int(abs(Rx))*iMagic
           iGy=int(abs(Ry))*iMagic
           iGz=int(abs(Rz))*iMagic
        endif
*
* make a stupid Patch: Cerius2 works well only with even nets!
*
      iGx=(iGx+1)/2*2
      iGy=(iGy+1)/2*2
      iGz=(iGz+1)/2*2
      mynCenter=nCenter

*----------------------------------------------------------------------*
*     Print coordinates of the system                                  *
*----------------------------------------------------------------------*
c      call bXML('Coord')
c      call iXML('nCoord',nCenter)
c      Do iAt=0,nCenter-1
c        Call cXML('Atom',AtomLbl(lw2+iAt))
c        call daXML('Atom coord',Work(ipCoor+3*iAt),3)
c      End Do
c      call eXML('Coord')

500   continue
      Write(6,*)
      Write(6,'(6X,A)')'Cartesian coordinates:'
      Write(6,'(6X,A)')'-----------------------------------------'
      Write(6,'(6X,A)')'No.  Label     X         Y         Z     '
      Write(6,'(6X,A)')'-----------------------------------------'
      Do iAt=0,nCenter-1
        Write(6,'(4X,I4,3X,A,2X,3F10.5)')
     &  iAt+1,AtomLbl(lw2+iAt),
     &  Work(ipCoor+3*iAt),Work(ipCoor+3*iAt+1),Work(ipCoor+3*iAt+2)
      End Do
      Write(6,'(6X,A)')'-----------------------------------------'
      Write(6,'(6X,A,3F12.6)')'Grid Origin      = ',Ox,Oy,Oz
      Write(6,'(6X,A,3F12.6)')'Grid Axis Length = ',Rx,Ry,Rz
      Write(6,'(6X,A)')'-----------------------------------------'
      Write(6,*)
      Write(6,*)

*
*----------------------------------------------------------------------*
*     Normal exit                                                      *
*----------------------------------------------------------------------*
      Return
      End
*----------------------------------------------------------------------*
