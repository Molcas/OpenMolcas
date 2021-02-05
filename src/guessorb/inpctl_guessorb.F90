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

use GuessOrb_Global, only: GapThr, iPrFmt, MxAtom, MxSym, PrintEor, PrintMOs, PrintPop, PrThr, SThr, TThr
use Definitions, only: iwp, u6

implicit none
!----------------------------------------------------------------------*
! Local data                                                           *
!----------------------------------------------------------------------*
logical(kind=iwp) :: Trace
character(len=180) :: Key, Line
integer(kind=iwp) :: LuSpool, itmp
!----------------------------------------------------------------------*
! External routines                                                    *
!----------------------------------------------------------------------*
integer(kind=iwp), external :: isFreeUnit
character(len=180), external :: Get_Ln
!----------------------------------------------------------------------*
! Setup                                                                *
!----------------------------------------------------------------------*
Trace = .false.
if (Trace) write(u6,*) '>>> Entering inpctl'
!----------------------------------------------------------------------*
! Process input                                                        *
!----------------------------------------------------------------------*
LuSpool = 17
LuSpool = isFreeUnit(LuSpool)
call SpoolInp(LuSpool)
call RdNLst(LuSpool,'GuessOrb')

input_loop: do
  Key = Get_Ln(LuSpool)
  Line = Key
  call UpCase(Line)
  select case (Line(1:4))

    !----------------------------------------------------------------------*
    ! NOMOs: skip printing of MOs, obsolete                                *
    !----------------------------------------------------------------------*
    case ('NOMO')
      write(u6,*) '******************************************'
      write(u6,*) '******************************************'
      write(u6,*) '***  OBSOLETE: do not use keyword NOMO ***'
      write(u6,*) '******************************************'
      write(u6,*) '******************************************'
      write(u6,*)
      PrintMOs = .false.

    !----------------------------------------------------------------------*
    ! PRMOs: MO print level.                                               *
    !----------------------------------------------------------------------*
    case ('PRMO')
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

    !----------------------------------------------------------------------*
    ! PRPOpulation: Mulliken print level                                   *
    !----------------------------------------------------------------------*
    case ('PRPO')
      PrintPop = .true.

    !----------------------------------------------------------------------*
    ! STHReshold: threshold for removing linear dependence, from S         *
    !----------------------------------------------------------------------*
    case ('STHR')
      Line = Get_Ln(LuSpool)
      call Get_F1(1,SThr)

    !----------------------------------------------------------------------*
    ! TTHReshold: threshold for removing linear dependence, from T         *
    !----------------------------------------------------------------------*
    case ('TTHR')
      Line = Get_Ln(LuSpool)
      call Get_F1(1,TThr)

    !----------------------------------------------------------------------*
    ! GapThr: threshold for homo-lumo gap.                                 *
    !----------------------------------------------------------------------*
    case ('GAPT')
      Line = Get_Ln(LuSpool)
      call Get_F1(1,GapThr)

    case ('END ')
      exit input_loop

    case default
      write(u6,*) 'InpCtl_GuessOrb: unidentified key word  : ',Key
      write(u6,*) 'InpCtl_GuessOrb: internal representation: ',Line(1:4)
      call FindErrorLine
      call Quit_OnUserError()
      exit input_loop
  end select
end do input_loop
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
if (Trace) write(u6,*) '<<< Exiting inpctl'

return

end subroutine InpCtl_GuessOrb
