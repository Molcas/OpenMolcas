************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2000, Markus P. Fuelscher                              *
************************************************************************
      Subroutine RdInp_FFPT
*
************************************************************************
*                                                                      *
*     Objective: Read and interprete input                             *
*                                                                      *
************************************************************************
*                                                                      *
*     In order to introduce the ffpt input to the CERIUS2 interface    *
*     the parser had to be changed!                                    *
*     M. Fuelscher, Lund Univeristy, Sweden, February 2000             *
*                                                                      *
************************************************************************
*
      Implicit Real*8 ( A-H,O-Z )
*

#include "input.fh"
*
      Parameter ( mCom = 11 )
      Character*20 FmtLog
      Character*4 Com(mCom)
      Data Com / 'TITL','DIPO','EFLD','QUAD','OCTU',
     &           'EFGR','RELA','GLBL','SELE','CUMU','END ' /
*
      Character*72 Line,Temp1,Temp2
      Character*4  Token
      Logical      Op0(9)
      Data         Op0 /9*.false./
      Logical      Op2(3)
      Data         Op2 /3*.false./
      Logical      Op3(4)
      Data         Op3 /4*.false./
      Logical      Op4(8)
      Data         Op4 /8*.false./
      Logical      Op5(12)
      Data         Op5 /12*.false./
      Logical      Op6(8)
      Data         Op6 /8*.false./

      LCumulate=.False.
*
*----------------------------------------------------------------------*
*                                                                      *
*     Start procedure                                                  *
*     Locate "start of input"                                          *
*                                                                      *
*----------------------------------------------------------------------*
*
      Call RdNlst(5,'FFPT')
      Temp2 = ' '
      Temp1 = ' '
      Line  = ' &FFPT &END'
*
*----------------------------------------------------------------------*
*     Initialize counters                                              *
*----------------------------------------------------------------------*
*
      newline = 0
      mLbl = 0
      mTit = 0
*
*----------------------------------------------------------------------*
*     Read the input stream line by line and identify key command      *
*----------------------------------------------------------------------*
*
      jCom=0
 1    Temp2 = Temp1
      Temp1 = Line
      Read(5,'(A)',Err=991,End=991) Line
      newline = newline+1
      Call LeftAd(Line)
      If ( Line(1:1).eq.' ' .or. Line(1:1).eq.'*' ) Goto 1
      Call StdFmt(Line,Token)
 2    jCom=0
      Do iCom=1,mCom
        If ( Token.eq.Com(iCom) ) jCom=iCom
      End Do
      If ( jCom.eq.0 ) Goto 992
*
*----------------------------------------------------------------------*
*     Branch to the processing of the command sections                 *
*----------------------------------------------------------------------*
*
      Goto (10,20,30,40,50,60,70,80,90,100,1000),jCom
*
*---  Process the "TITL" command --------------------------------------*
*
 10   Continue
      If ( Op0(1) ) Goto 993
      Op0(1) = .true.
      ComStk(1,0,0,0)=.true.
 15   Temp2 = Temp1
      Temp1 = Line
      newline = newline+1
      Read(5,'(A)',Err=991,End=991) Line
      Call LeftAd(Line)
      If ( Line(1:1).eq.' ' .or. Line(1:1).eq.'*' ) Goto 15
      Call StdFmt(Line,Token)
      Do iCom=1,mCom
        If ( Token.eq.Com(iCom) ) Goto 2
      End Do
      mTit=mTit+1
      If ( mTit.le.MxTitL ) then
         Title(mTit)=Line
      End If
      Goto 15
*
*---  Process the "DIPO" command --------------------------------------*
*
 20   Continue
      If ( Op0(2) ) Goto 993
      Op0(2) = .true.
 25   Temp2 = Temp1
      Temp1 = Line
      newline = newline+1
      Read(5,'(A)',Err=991,End=991) Line
      Call LeftAd(Line)
      If ( Line(1:1).eq.' ' .or. Line(1:1).eq.'*' ) Goto 25
      Call UpCase(Line)
      Do i = 1,LEN(Line)
        If ( Line(i:i).eq.'=' ) Line(i:i) = ' '
        If ( Line(i:i).eq.CHAR(9) ) Line(i:i) = ' '
      End Do
      Call StdFmt(Line,Token)
      Do iCom=1,mCom
        If ( Token.eq.Com(iCom) ) Goto 2
      End Do
      i1 = 2
      i2 = LEN(Line)
      If ( Token.eq.'X' ) then
        If ( Op2(1) ) Goto 993
        Op2(1) = .true.
        Read(Line(i1:i2),*,Err=991,End=991) W
        ComStk(2,1,1,0) = .true.
        ComStk(2,1,1,1) = .true.
        ComVal(2,1,1,1) = W
        Goto 25
      End If
      If ( Token.eq.'Y' ) then
        If ( Op2(2) ) Goto 993
        Op2(2) = .true.
        Read(Line(i1:i2),*,Err=991,End=991) W
        ComStk(2,1,1,0) = .true.
        ComStk(2,1,1,2) = .true.
        ComVal(2,1,1,2) = W
        Goto 25
      End If
      If ( Token.eq.'Z' ) then
        If ( Op2(3) ) Goto 993
        Op2(3) = .true.
        Read(Line(i1:i2),*,Err=991,End=991) W
        ComStk(2,1,1,0) = .true.
        ComStk(2,1,1,3) = .true.
        ComVal(2,1,1,3) = W
        Goto 25
      End If
      Goto 992
*
*---  Process the "EFLD" command --------------------------------------*
*
 30   Continue
      If ( Op0(3) ) Goto 993
      Op0(3) = .true.
 35   Temp2 = Temp1
      Temp1 = Line
      newline = newline+1
      Read(5,'(A)',Err=991,End=991) Line
      Call LeftAd(Line)
      If ( Line(1:1).eq.' ' .or. Line(1:1).eq.'*' ) Goto 35
      Call UpCase(Line)
      Do i = 1,LEN(Line)
        If ( Line(i:i).eq.'=' ) Line(i:i) = ' '
        If ( Line(i:i).eq.CHAR(9) ) Line(i:i) = ' '
      End Do
      Call StdFmt(Line,Token)
      Do iCom=1,mCom
        If ( Token.eq.Com(iCom) ) Goto 2
      End Do
      i1 = 2
      i2 = LEN(Line)
      If ( Token.eq.'X' ) then
        If ( Op3(1) ) Goto 993
        Op3(1) = .true.
        Read(Line(i1:i2),*,Err=991,End=991) W
        ComStk(2,3,1,0) = .true.
        ComStk(2,3,1,1) = .true.
        ComVal(2,3,1,1) = W
        Goto 35
      End If
      If ( Token.eq.'Y' ) then
        If ( Op3(2) ) Goto 993
        Op3(2) = .true.
        Read(Line(i1:i2),*,Err=991,End=991) W
        ComStk(2,3,1,0) = .true.
        ComStk(2,3,1,2) = .true.
        ComVal(2,3,1,2) = W
        Goto 35
      End If
      If ( Token.eq.'Z' ) then
        If ( Op3(3) ) Goto 993
        Op3(3) = .true.
        Read(Line(i1:i2),*,Err=991,End=991) W
        ComStk(2,3,1,0) = .true.
        ComStk(2,3,1,3) = .true.
        ComVal(2,3,1,3) = W
        Goto 35
      End If
      i1 = 5
      If ( Token.eq.'ORIG' ) then
        If ( Op3(4) ) Goto 993
        Op3(4) = .true.
        i1=5
        Read(Line(i1:i2),*,Err=991,End=991) X,Y,Z
        ComStk(2,3,2,0) = .true.
        ComStk(2,3,2,1) = .true.
        ComVal(2,3,2,1) = X
        ComStk(2,3,2,2) = .true.
        ComVal(2,3,2,2) = Y
        ComStk(2,3,2,3) = .true.
        ComVal(2,3,2,3) = Z
        Goto 35
      End If
      Goto 992
*
*---  Process the "QUAD" command --------------------------------------*
*
 40   Continue
      If ( Op0(4) ) Goto 993
      Op0(4) = .true.
 45   Temp2 = Temp1
      Temp1 = Line
      newline = newline+1
      Read(5,'(A)',End=991) Line
      Call LeftAd(Line)
      If ( Line(1:1).eq.' ' .or. Line(1:1).eq.'*' ) Goto 45
      Call UpCase(Line)
      Do i = 1,LEN(Line)
        If ( Line(i:i).eq.'=' ) Line(i:i) = ' '
        If ( Line(i:i).eq.CHAR(9) ) Line(i:i) = ' '
      End Do
      Call StdFmt(Line,Token)
      Do iCom=1,mCom
        If ( Token.eq.Com(iCom) ) Goto 2
      End Do
      i1 = 3
      i2 = LEN(Line)
      If ( Token.eq.'XX' ) then
        If ( Op4(1) ) Goto 993
        Op4(1) = .true.
        Read(Line(i1:i2),*,Err=991,End=991) W
        ComStk(2,2,1,0) = .true.
        ComStk(2,2,1,1) = .true.
        ComVal(2,2,1,1) = W
        Goto 45
      End If
      If ( Token.eq.'YY' ) then
        If ( Op4(2) ) Goto 993
        Op4(2) = .true.
        Read(Line(i1:i2),*,Err=991,End=991) W
        ComStk(2,2,1,0) = .true.
        ComStk(2,2,1,4) = .true.
        ComVal(2,2,1,4) = W
        Goto 45
      End If
      If ( Token.eq.'ZZ' ) then
        If ( Op4(3) ) Goto 993
        Op4(3) = .true.
        Read(Line(i1:i2),*,Err=991,End=991) W
        ComStk(2,2,1,0) = .true.
        ComStk(2,2,1,6) = .true.
        ComVal(2,2,1,6) = W
        Goto 45
      End If
      If ( Token.eq.'XY' ) then
        If ( Op4(4) ) Goto 993
        Op4(4) = .true.
        Read(Line(i1:i2),*,Err=991,End=991) W
        ComStk(2,2,1,0) = .true.
        ComStk(2,2,1,2) = .true.
        ComVal(2,2,1,2) = W
        Goto 45
      End If
      If ( Token.eq.'XZ' ) then
        If ( Op4(5) ) Goto 993
        Op4(5) = .true.
        Read(Line(i1:i2),*,Err=991,End=991) W
        ComStk(2,2,1,0) = .true.
        ComStk(2,2,1,3) = .true.
        ComVal(2,2,1,3) = W
        Goto 45
      End If
      If ( Token.eq.'YZ' ) then
        If ( Op4(6) ) Goto 993
        Op4(6) = .true.
        Read(Line(i1:i2),*,Err=991,End=991) W
        ComStk(2,2,1,0) = .true.
        ComStk(2,2,1,5) = .true.
        ComVal(2,2,1,5) = W
        Goto 45
      End If
      If ( Token.eq.'RR' ) then
        If ( Op4(7) ) Goto 993
        Op4(7) = .true.
        Read(Line(i1:i2),*,Err=991,End=991) W
        ComStk(2,2,1,0) = .true.
        ComStk(2,2,1,7) = .true.
        ComVal(2,2,1,7) = W
        Goto 45
      End If
      If ( Token.eq.'ORIG' ) then
        If ( Op4(8) ) Goto 993
        Op4(8) = .true.
        i1=5
        Read(Line(i1:i2),*,Err=991,End=991) X,Y,Z
        ComStk(2,2,2,0) = .true.
        ComStk(2,2,2,1) = .true.
        ComVal(2,2,2,1) = X
        ComStk(2,2,2,2) = .true.
        ComVal(2,2,2,2) = Y
        ComStk(2,2,2,3) = .true.
        ComVal(2,2,2,3) = Z
        Goto 45
      End If
      Goto 992
*
*---  Process the "OCTU" command --------------------------------------*
*
 50   Continue
      If ( Op0(5) ) Goto 993
      Op0(5) = .true.
 55   Temp2 = Temp1
      Temp1 = Line
      newline = newline+1
      Read(5,'(A)',Err=991,End=991) Line
      Call LeftAd(Line)
      If ( Line(1:1).eq.' ' .or. Line(1:1).eq.'*' ) Goto 55
      Call UpCase(Line)
      Do i = 1,LEN(Line)
        If ( Line(i:i).eq.'=' ) Line(i:i) = ' '
        If ( Line(i:i).eq.CHAR(9) ) Line(i:i) = ' '
      End Do
      Call StdFmt(Line,Token)
      Do iCom=1,mCom
        If ( Token.eq.Com(iCom) ) Goto 2
      End Do
      i1 = 4
      i2 = LEN(Line)
      If ( Token.eq.'XXX' ) then
        If ( Op5(1) ) Goto 993
        Op5(1) = .true.
        Read(Line(i1:i2),*,Err=991,End=991) W
        ComStk(2,6,1,1) = .true.
        ComStk(2,6,1,1) = .true.
        ComVal(2,6,1,1) = W
        Goto 55
      End If
      If ( Token.eq.'XYY' ) then
        If ( Op5(2) ) Goto 993
        Op5(2) = .true.
        Read(Line(i1:i2),*,Err=991,End=991) W
        ComStk(2,6,1,4) = .true.
        ComStk(2,6,1,4) = .true.
        ComVal(2,6,1,4) = W
        Goto 55
      End If
      If ( Token.eq.'XZZ' ) then
        If ( Op5(3) ) Goto 993
        Op5(3) = .true.
        Read(Line(i1:i2),*,Err=991,End=991) W
        ComStk(2,6,1,6) = .true.
        ComStk(2,6,1,6) = .true.
        ComVal(2,6,1,6) = W
        Goto 55
      End If
      If ( Token.eq.'XXY' ) then
        If ( Op5(4) ) Goto 993
        Op5(4) = .true.
        Read(Line(i1:i2),*,Err=991,End=991) W
        ComStk(2,6,1,2) = .true.
        ComStk(2,6,1,2) = .true.
        ComVal(2,6,1,2) = W
        Goto 55
      End If
      If ( Token.eq.'YYY' ) then
        If ( Op5(5) ) Goto 993
        Op5(5) = .true.
        Read(Line(i1:i2),*,Err=991,End=991) W
        ComStk(2,6,1,7) = .true.
        ComStk(2,6,1,7) = .true.
        ComVal(2,6,1,7) = W
        Goto 55
      End If
      If ( Token.eq.'YZZ' ) then
        If ( Op5(6) ) Goto 993
        Op5(6) = .true.
        Read(Line(i1:i2),*,Err=991,End=991) W
        ComStk(2,6,1,9) = .true.
        ComStk(2,6,1,9) = .true.
        ComVal(2,6,1,9) = W
        Goto 55
      End If
      If ( Token.eq.'XXZ' ) then
        If ( Op5(7) ) Goto 993
        Op5(8) = .true.
        Read(Line(i1:i2),*,Err=991,End=991) W
        ComStk(2,6,1,3) = .true.
        ComStk(2,6,1,3) = .true.
        ComVal(2,6,1,3) = W
        Goto 55
      End If
      If ( Token.eq.'YYZ' ) then
        If ( Op5(9) ) Goto 993
        Op5(9) = .true.
        Read(Line(i1:i2),*,Err=991,End=991) W
        ComStk(2,6,1,8) = .true.
        ComStk(2,6,1,8) = .true.
        ComVal(2,6,1,8) = W
        Goto 55
      End If
      If ( Token.eq.'ZZZ' ) then
        If ( Op5(10) ) Goto 993
        Op5(10) = .true.
        Read(Line(i1:i2),*,Err=991,End=991) W
        ComStk(2,6,1,10) = .true.
        ComStk(2,6,1,10) = .true.
        ComVal(2,6,1,10) = W
        Goto 55
      End If
      If ( Token.eq.'XYZ' ) then
        If ( Op5(11) ) Goto 993
        Op5(11) = .true.
        Read(Line(i1:i2),*,Err=991,End=991) W
        ComStk(2,6,1,5) = .true.
        ComStk(2,6,1,5) = .true.
        ComVal(2,6,1,5) = W
        Goto 55
      End If
      If ( Token.eq.'ORIG' ) then
        If ( Op5(12) ) Goto 993
        Op5(12) = .true.
        i1=5
        Read(Line(i1:i2),*,Err=991,End=991) X,Y,Z
        ComStk(2,6,2,0) = .true.
        ComStk(2,6,2,1) = .true.
        ComVal(2,6,2,1) = X
        ComStk(2,6,2,2) = .true.
        ComVal(2,6,2,2) = Y
        ComStk(2,6,2,3) = .true.
        ComVal(2,6,2,3) = Z
        Goto 55
      End If
      Goto 992
*
*---  Process the "EFGR" command --------------------------------------*
*
 60   Continue
      If ( Op0(6) ) Goto 993
      Op0(6) = .true.
 65   Temp2 = Temp1
      Temp1 = Line
      newline = newline+1
      Read(5,'(A)',Err=991,End=991) Line
      Call LeftAd(Line)
      If ( Line(1:1).eq.' ' .or. Line(1:1).eq.'*' ) Goto 65
      Call UpCase(Line)
      Do i = 1,LEN(Line)
        If ( Line(i:i).eq.'=' ) Line(i:i) = ' '
        If ( Line(i:i).eq.CHAR(9) ) Line(i:i) = ' '
      End Do
      Call StdFmt(Line,Token)
      Do iCom=1,mCom
        If ( Token.eq.Com(iCom) ) Goto 2
      End Do
      i1 = 3
      i2 = LEN(Line)
      If ( Token.eq.'XX' ) then
        If ( Op6(1) ) Goto 993
        Op6(1) = .true.
        Read(Line(i1:i2),*,Err=991,End=991) W
        ComStk(2,4,1,0) = .true.
        ComStk(2,4,1,1) = .true.
        ComVal(2,4,1,1) = W
        Goto 65
      End If
      If ( Token.eq.'YY' ) then
        If ( Op6(2) ) Goto 993
        Op6(2) = .true.
        Read(Line(i1:i2),*,Err=991,End=991) W
        ComStk(2,4,1,0) = .true.
        ComStk(2,4,1,4) = .true.
        ComVal(2,4,1,4) = W
        Goto 65
      End If
      If ( Token.eq.'ZZ' ) then
        If ( Op6(3) ) Goto 993
        Op6(3) = .true.
        Read(Line(i1:i2),*,Err=991,End=991) W
        ComStk(2,4,1,0) = .true.
        ComStk(2,4,1,6) = .true.
        ComVal(2,4,1,6) = W
        Goto 65
      End If
      If ( Token.eq.'XY' ) then
        If ( Op6(4) ) Goto 993
        Op6(4) = .true.
        Read(Line(i1:i2),*,Err=991,End=991) W
        ComStk(2,4,1,0) = .true.
        ComStk(2,4,1,2) = .true.
        ComVal(2,4,1,2) = W
        Goto 65
      End If
      If ( Token.eq.'XZ' ) then
        If ( Op6(5) ) Goto 993
        Op6(5) = .true.
        Read(Line(i1:i2),*,Err=991,End=991) W
        ComStk(2,4,1,0) = .true.
        ComStk(2,4,1,3) = .true.
        ComVal(2,4,1,3) = W
        Goto 65
      End If
      If ( Token.eq.'YZ' ) then
        If ( Op6(6) ) Goto 993
        Op6(6) = .true.
        Read(Line(i1:i2),*,Err=991,End=991) W
        ComStk(2,4,1,0) = .true.
        ComStk(2,4,1,5) = .true.
        ComVal(2,4,1,5) = W
        Goto 65
      End If
      If ( Token.eq.'ORIG' ) then
        If ( Op6(8) ) Goto 993
        Op6(8) = .true.
        i1=5
        Read(Line(i1:i2),*,Err=991,End=991) X,Y,Z
        ComStk(2,4,2,0) = .true.
        ComStk(2,4,2,1) = .true.
        ComVal(2,4,2,1) = X
        ComStk(2,4,2,2) = .true.
        ComVal(2,4,2,2) = Y
        ComStk(2,4,2,3) = .true.
        ComVal(2,4,2,3) = Z
        Goto 65
      End If
      Goto 992
*
*---  Process the "RELA" command --------------------------------------*
*
 70   Continue
      If ( Op0(7) ) Goto 993
      Op0(7) = .true.
 75   Temp2 = Temp1
      Temp1 = Line
      newline = newline+1
      Read(5,'(A)',Err=991,End=991) Line
      Call LeftAd(Line)
      If ( Line(1:1).eq.' ' .or. Line(1:1).eq.'*' ) Goto 75
      i1 = 1
      i2 = LEN(Line)
      Read(Line(i1:i2),*,Err=991,End=991) W
      ComStk(2,5,0,0) = .true.
      ComStk(2,5,0,1) = .true.
      ComVal(2,5,0,1) = W
      Goto 1
*
*---  Process the "GLBL" command --------------------------------------*
*
 80   Continue
      If ( Op0(8) ) Goto 993
      Op0(8) = .true.
      ComStk(3,0,0,0)=.true.
 85   Temp2 = Temp1
      Temp1 = Line
      newline = newline+1
      Read(5,'(A)',Err=991,End=991) Line
      Call LeftAd(Line)
      If ( Line(1:1).eq.' ' .or. Line(1:1).eq.'*' ) Goto 85
      Call UpCase(Line)
      Call StdFmt(Line,Token)
      Do iCom=1,mCom
        If ( Token.eq.Com(iCom) ) Goto 2
      End Do
      i1 = 1
      i2 = LEN(Line)
      mLbl = mLbl+1
      If ( mLbl.gt.MxLbl ) Goto 994
      iSta = INDEX(Line,'''')
      Line(i1:iSta) = ' '
      iEnd = INDEX(Line,'''')
      gLblN(mLbl) = Line(iSta+1:iEnd-1)
      Line(i1:iEnd) = ' '
      Read(Line(i1:i2),*,End=991) C,W
      gLblC(mLbl) = NINT(C)
      gLblW(mLbl) = W
      Do i = 1,mLbl-1
        If ( gLblN(i).eq.gLblN(mLbl)
     &      .and.
     &       gLblC(i).eq.gLblC(mLbl) ) Goto 993
      End Do
      Goto 85
*
*---  Process the "SELE" command --------------------------------------*
*
 90   Continue
*-- Initialize
      Do i=1,MxSets
        Atoms(i)=.false.
        Do j=1,MxSets
          Bonds(i,j)=.false.
        Enddo
      Enddo
      ComStk(4,0,0,0)=.true.
 95   Temp2=Temp1
      Temp1=Line
      newline=newline+1
      Read(5,'(A)',Err=991,End=991) Line
      Call LeftAd(Line)
      If ( Line(1:1).eq.' ' .or. Line(1:1).eq.'*' ) Goto 95
      Read(Line,*)nSets
      Do i=1,nSets
        Read(5,*)Atoms(i),iSelection(1,i),iSelection(2,i)
      Enddo
      Do i=2,nSets
        If(i.lt.10) Write(FmtLog,79121)'(',i-1,'L2)'
        If(i.ge.10) Write(FmtLog,79122)'(',i-1,'L2)'
        Read(5,FmtLog)(Bonds(i,j),j=1,i-1)
        Do j=1,i-1
          Bonds(j,i)=Bonds(i,j)
        Enddo
      Enddo
      Read(5,*)(TranCoo(k),k=1,3)
*      Read(5,*)SiffBond
79121 Format(A,I1,A)
79122 Format(A,I2,A)
      Go to 1
*
*---  Process the "CUMU" command --------------------------------------*
*
 100  Continue
*     Add the perturbation to the current H0 instead of
*     to the vacuum H0. This enables multiple FFPT runs after eachother,
*     useful when using SELE.
      LCumulate=.True.
      Go to 1

*
*---  Process the "END " command --------------------------------------*
*
 1000 Continue
      If ( Op0(9) ) Goto 993
      Op0(9) = .true.
*      ComStk(4,0,0,0) = .true.
      ComStk(5,0,0,0) = .true.
*---  Check for redunancy in the origin input
      If ( ComStk(2,2,2,4) ) Then
         If ( ComStk(2,2,2,1) .or.
     &        ComStk(2,2,2,2) .or.
     &        ComStk(2,2,2,3)      ) Goto 996
      End If
      If ( ComStk(2,3,2,4) ) Then
         If ( ComStk(2,3,2,1) .or.
     &        ComStk(2,3,2,2) .or.
     &        ComStk(2,3,2,3)      ) Goto 996
      End If
      If ( ComStk(2,4,2,4) ) Then
         If ( ComStk(2,4,2,1) .or.
     &        ComStk(2,4,2,2) .or.
     &        ComStk(2,4,2,3)      ) Goto 996
      End If
      If ( ComStk(2,6,2,4) ) Then
         If ( ComStk(2,6,2,1) .or.
     &        ComStk(2,6,2,2) .or.
     &        ComStk(2,6,2,3)      ) Goto 996
      End If
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
*
*----------------------------------------------------------------------*
*     Error handling                                                   *
*----------------------------------------------------------------------*
*
 991  Write (6,*)
      Write(6,'(2X,A)') 'The program failed to read the input.'
      Write(6,'(2X,A)') 'Please check your input data.'
      Write(6,*)
      Write(6,'(2X,A,I3.3,A)') 'The error occured at line',newline,
     &                         ' after the &FFPT &END line'
      Write(6,'(2X,A,A)') 'The current line is:      ',Line
      Write(6,'(2X,A,A)') 'The previous line is:     ',Temp1
      Write(6,'(2X,A,A)') 'The next previous line is:',Temp2
      Call Quit_OnUserError
 992  Write (6,*)
      Write(6,'(2X,A)') 'The program has been supplied with an unknown'
      Write(6,'(2X,A)') 'keyword. Please correct your input data.'
      Write(6,*)
      Write(6,'(2X,A,I3.3,A)') 'The error occured at line',newline,
     &                         ' after the &FFPT &END line'
      Write(6,'(2X,A,A)') 'The current line is:      ',Line
      Write(6,'(2X,A,A)') 'The previous line is:     ',Temp1
      Write(6,'(2X,A,A)') 'The next previous line is:',Temp2
      Call Quit_OnUserError
 993  Write (6,*)
      Write(6,'(2X,A)') 'A command or one of its components has been'
      Write(6,'(2X,A)') 'multiply defined. Please correct your input.'
      Write(6,*)
      Write(6,'(2X,A,I3.3,A)') 'The error occured at line',newline,
     &                         ' after the &FFPT &END line'
      Write(6,'(2X,A,A)') 'The current line is:      ',Line
      Write(6,'(2X,A,A)') 'The previous line is:     ',Temp1
      Write(6,'(2X,A,A)') 'The next previous line is:',Temp2
      Call Quit_OnUserError
 994  Write (6,*)
      Write(6,'(2X,A)') 'The number of perturbations requested exceeds'
      Write(6,'(2X,A)') 'the internal buffer size. Increase the para-'
      Write(6,'(2X,A)') 'meter MxLbl and recompile the program.'
      Call Quit_OnUserError
 996  Write (6,*)
      Write(6,'(2X,A)') 'The definition of the origin of an operator '
      Write(6,'(2X,A)') 'is not unique.  Please correct your input.'
      Call Quit_OnUserError
      End
