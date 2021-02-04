!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2004, Per-Olof Widmark                                 *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This routine read the input for guessorb.                            *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: Oct 2004                                                    *
!                                                                      *
!***********************************************************************
      Subroutine InpCtl_GuessOrb()
use GuessOrb_Global, only: GapThr, iPrFmt, LenIn, LenIn1, LenIn8, MxAtom, MxSym, PrintEor, PrintMOs, PrintPop, PrThr, SThr, TThr
      Implicit None
!----------------------------------------------------------------------*
! Parameters                                                           *
!----------------------------------------------------------------------*
      Character*15 myName
      Parameter (myName='InpCtl_GuessOrb')
!----------------------------------------------------------------------*
! Local data                                                           *
!----------------------------------------------------------------------*
      Logical       Trace
      Character*180 Key, Line
      Character*180 Get_Ln
      External      Get_Ln
      Integer       LuSpool
      Integer       isFreeUnit
      Integer       itmp
!----------------------------------------------------------------------*
! External routines                                                    *
!----------------------------------------------------------------------*
      External      isFreeUnit
!----------------------------------------------------------------------*
! Setup                                                                *
!----------------------------------------------------------------------*
      Trace=.false.
      If(Trace) Write(6,*) '>>> Entering inpctl'
!----------------------------------------------------------------------*
! Process input                                                        *
!----------------------------------------------------------------------*
      LuSpool=17
      LuSpool=isFreeUnit(LuSpool)
      Call SpoolInp(LuSpool)
      Call RdNLst(LuSpool,'GuessOrb')

  999 Continue
      Key=Get_Ln(LuSpool)
      Line=Key
      Call UpCase(Line)
      If (Line(1:4).eq.'NOMO') Go To  1000
      If (Line(1:4).eq.'PRMO') Go To  1100
      If (Line(1:4).eq.'PRPO') Go To  1200
      If (Line(1:4).eq.'STHR') Go To  1300
      If (Line(1:4).eq.'TTHR') Go To  1400
      If (Line(1:4).eq.'GAPT') Go To  1500
      If (Line(1:4).eq.'END ') Go To 99999
      Write(6,*) myName,': unidentified key word  : ',Key
      Write(6,*) myName,': internal representation: ',Line(1:4)
      Call FindErrorLine
      Call Quit_OnUserError()
!----------------------------------------------------------------------*
! NOMOs: skip printing of MOs, obsolete                                *
!----------------------------------------------------------------------*
 1000 Continue
      Write(6,*) '******************************************'
      Write(6,*) '******************************************'
      Write(6,*) '***  OBSOLETE: do not use keyword NOMO ***'
      Write(6,*) '******************************************'
      Write(6,*) '******************************************'
      Write(6,*)
      PrintMOs=.False.
      Go To 999
!----------------------------------------------------------------------*
! PRMOs: MO print level.                                               *
!----------------------------------------------------------------------*
 1100 Continue
      Line=Get_Ln(LuSpool)
      Line(178:180)='5.0'
      Call Put_Ln(Line)
      Call Get_I1(1,itmp)
      Call Get_F1(2,PrThr)
      If(itmp.ge.4) Then
         PrintMOs=.true.
         PrintEor=.true.
         iPrFmt=3
      Else If(itmp.eq.3) Then
         PrintMOs=.true.
         PrintEor=.true.
         iPrFmt=2
      Else If(itmp.eq.2) Then
         PrintMOs=.true.
         PrintEor=.true.
         iPrFmt=1
      Else If(itmp.eq.1) Then
         PrintMOs=.true.
         PrintEor=.false.
         iPrFmt=1
      Else
         PrintMOs=.False.
         PrintEor=.false.
      End If
      Go To 999
!----------------------------------------------------------------------*
! PRPOpulation: Mulliken print level                                   *
!----------------------------------------------------------------------*
 1200 Continue
      PrintPop=.true.
      Go To 999
!----------------------------------------------------------------------*
! STHReshold: threshold for removing linear dependence, from S         *
!----------------------------------------------------------------------*
 1300 Continue
      Line=Get_Ln(LuSpool)
      Call Get_F1(1,SThr)
      Go To 999
!----------------------------------------------------------------------*
! TTHReshold: threshold for removing linear dependence, from T         *
!----------------------------------------------------------------------*
 1400 Continue
      Line=Get_Ln(LuSpool)
      Call Get_F1(1,TThr)
      Go To 999
!----------------------------------------------------------------------*
! GapThr: threshold for homo-lumo gap.                                 *
!----------------------------------------------------------------------*
 1500 Continue
      Line=Get_Ln(LuSpool)
      Call Get_F1(1,GapThr)
      Go To 999
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
99999 Continue
      If(Trace) Write(6,*) '<<< Exiting inpctl'
      Return
      End
