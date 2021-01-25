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
      Subroutine BasisReader(LuWr,nBase,
     &      iglobal,nxbas,xb_label,xb_bas,iErr)
      Implicit Integer (i-n)
      Implicit Real*8 (a-h,o-z)
      Character*48 Line
      Character*2  SimbA
#include "g_zmatconv.fh"
      Logical Found
        character *(*) xb_label(*)
        character *(*) xb_bas(*)
      character *48 TMP
      iErr = 0
      nBase = 0
      icount=1

10    CONTINUE
C Read(iZMUnit,'(A)',Err=9904,End=9999) Line
      if(iglobal.eq.0) then
        Line=trim(xb_label(icount))//'.'//trim(xb_bas(icount))
      else
              Line='FF.'//trim(xb_bas(icount))
        nxbas=1
      endif
c      print *,'>>',Line
      Found = .False.
      SimbA = Line(1:2)
      If (SimbA(2:2).EQ.'.') SimbA(2:2)=' '
      Do i = 1, Num_Elem
        If ((SimbA(1:1).le.'z').and.(SimbA(1:1).ge.'a'))
     &       SimbA(1:1)=char(ichar(SimbA(1:1))-32)
        If ((SimbA(2:2).le.'Z').and.(SimbA(2:2).ge.'A'))
     &       SimbA(2:2)=char(ichar(SimbA(2:2))+32)

      If (SimbA.EQ.AdjustL(PTab(i)).and.iglobal.eq.0) then
          Base(i) = Line
          BasAva(i) = .True.
          Found = .True.
          nBase = nBase + 1
      EndIf
      If (iglobal.eq.1) then
          Base(i) = Line
          Base(i)(1:2) = AdjustL(PTab(i))
          if(Base(i)(2:2).eq.' ') then
            TMP=Base(i)(3:)
            Base(i)(2:)=TMP(1:47)
          endif
          BasAva(i) = .True.
          Found = .True.
          nBase = nBase + 1
      EndIf

      EndDo
      If (.NOT.Found) then
        iErr = 1
        Write(LuWr,*) ' [BasisReader]: Wrong symbol in line'
        Write(LuWr,*) '                ',Line
        GoTo 9999
      EndIf
      icount=icount+1
      if(icount.le.nxbas) GoTo 10


c9904  iErr = 1
c      Write(LuWr,*) ' [BasisReader]: Unable to read z-matrix !'
c      GoTo 9999
9999  continue
      Return
      End


      Subroutine BasisConsistency(LuWr,iErr)
      Implicit Integer (i-n)
      Implicit Real*8 (a-h,o-z)
#include "g_zmatconv.fh"

      iErr = 0
      Do i = 1, 100
        If (BasReq(i) .AND. .NOT.BasAva(i)) GoTo 9900
      EndDo
      Return

9900  iErr = 1
      Write(LuWr,*) ' [BasisConsistency]: Atom NA=',i,' requires BS'
      Return
      End
