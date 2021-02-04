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

subroutine InpCtl_GuessOrb()

use GuessOrb_Global, only: GapThr, iPrFmt, LenIn, LenIn1, LenIn8, MxAtom, MxSym, PrintEor, PrintMOs, PrintPop, PrThr, SThr, TThr

implicit none
!----------------------------------------------------------------------*
! Parameters                                                           *
!----------------------------------------------------------------------*
character*15 myName
parameter(myName='InpCtl_GuessOrb')
!----------------------------------------------------------------------*
! Local data                                                           *
!----------------------------------------------------------------------*
logical Trace
character*180 Key, Line
character*180 Get_Ln
external Get_Ln
integer LuSpool
integer isFreeUnit
integer itmp
!----------------------------------------------------------------------*
! External routines                                                    *
!----------------------------------------------------------------------*
external isFreeUnit
!----------------------------------------------------------------------*
! Setup                                                                *
!----------------------------------------------------------------------*
Trace = .false.
if (Trace) write(6,*) '>>> Entering inpctl'
!----------------------------------------------------------------------*
! Process input                                                        *
!----------------------------------------------------------------------*
LuSpool = 17
LuSpool = isFreeUnit(LuSpool)
call SpoolInp(LuSpool)
call RdNLst(LuSpool,'GuessOrb')

999 continue
Key = Get_Ln(LuSpool)
Line = Key
call UpCase(Line)
if (Line(1:4) == 'NOMO') Go To 1000
if (Line(1:4) == 'PRMO') Go To 1100
if (Line(1:4) == 'PRPO') Go To 1200
if (Line(1:4) == 'STHR') Go To 1300
if (Line(1:4) == 'TTHR') Go To 1400
if (Line(1:4) == 'GAPT') Go To 1500
if (Line(1:4) == 'END ') Go To 99999
write(6,*) myName,': unidentified key word  : ',Key
write(6,*) myName,': internal representation: ',Line(1:4)
call FindErrorLine
call Quit_OnUserError()
!----------------------------------------------------------------------*
! NOMOs: skip printing of MOs, obsolete                                *
!----------------------------------------------------------------------*
1000 continue
write(6,*) '******************************************'
write(6,*) '******************************************'
write(6,*) '***  OBSOLETE: do not use keyword NOMO ***'
write(6,*) '******************************************'
write(6,*) '******************************************'
write(6,*)
PrintMOs = .false.
Go To 999
!----------------------------------------------------------------------*
! PRMOs: MO print level.                                               *
!----------------------------------------------------------------------*
1100 continue
Line = Get_Ln(LuSpool)
Line(178:180) = '5.0'
call Put_Ln(Line)
call Get_I1(1,itmp)
call Get_F1(2,PrThr)
if (itmp >= 4) then
  PrintMOs = .true.
  PrintEor = .true.
  iPrFmt = 3
else if (itmp == 3) then
  PrintMOs = .true.
  PrintEor = .true.
  iPrFmt = 2
else if (itmp == 2) then
  PrintMOs = .true.
  PrintEor = .true.
  iPrFmt = 1
else if (itmp == 1) then
  PrintMOs = .true.
  PrintEor = .false.
  iPrFmt = 1
else
  PrintMOs = .false.
  PrintEor = .false.
end if
Go To 999
!----------------------------------------------------------------------*
! PRPOpulation: Mulliken print level                                   *
!----------------------------------------------------------------------*
1200 continue
PrintPop = .true.
Go To 999
!----------------------------------------------------------------------*
! STHReshold: threshold for removing linear dependence, from S         *
!----------------------------------------------------------------------*
1300 continue
Line = Get_Ln(LuSpool)
call Get_F1(1,SThr)
Go To 999
!----------------------------------------------------------------------*
! TTHReshold: threshold for removing linear dependence, from T         *
!----------------------------------------------------------------------*
1400 continue
Line = Get_Ln(LuSpool)
call Get_F1(1,TThr)
Go To 999
!----------------------------------------------------------------------*
! GapThr: threshold for homo-lumo gap.                                 *
!----------------------------------------------------------------------*
1500 continue
Line = Get_Ln(LuSpool)
call Get_F1(1,GapThr)
Go To 999
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
99999 continue
if (Trace) write(6,*) '<<< Exiting inpctl'

return

end subroutine InpCtl_GuessOrb
