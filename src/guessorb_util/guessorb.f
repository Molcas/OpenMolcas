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
* Copyright (C) 2004, Per-Olof Widmark                                 *
************************************************************************
      subroutine guessorb(iReturn,StandAlone)
************************************************************************
*                                                                      *
* This program computes start orbitals for use in SCF/RASSCF etc.      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University                                             *
*          Sweden                                                      *
* Written: Oct 2004                                                    *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Implementation notes:                                                *
*                                                                      *
************************************************************************
      Implicit None
#include "Molcas.fh"
#include "commgo.fh"
*----------------------------------------------------------------------*
* Dummy arguments                                                      *
*----------------------------------------------------------------------*
      Integer iReturn
      Logical StandAlone
*----------------------------------------------------------------------*
* Local variables.                                                     *
*----------------------------------------------------------------------*
      Integer iRC, iUHF
*----------------------------------------------------------------------*
* Prologue                                                             *
*----------------------------------------------------------------------*
      iReturn=0
      Call qEnter('guessorb')
      Call InitGO(StandAlone)
      If(StandAlone) Call InpCtl_GuessOrb
*----------------------------------------------------------------------*
* Select method to be used.                                            *
*----------------------------------------------------------------------*
      Call cre_gsswfn
      Call FckByInt(iRC,StandAlone)
*     If(iRC.eq.0) GoTo 999
      If(.true.) GoTo 999
      If(nSym.eq.1) Then
         Call Fmod1n(StandAlone)
      Else
         Call Fmod1s(StandAlone)
      End If
999   Continue
      Call cls_gsswfn
*----------------------------------------------------------------------*
* Produce MOLDEN input                                                 *
*----------------------------------------------------------------------*
      iUHF=0
      If (iRC.eq.0) then
        Call Molden_Interface(iUHF,'GSSORB','MD_GSS')
c        call grid_driver(-1,'SEWARD','GSSORB',iRc)
      End If
*----------------------------------------------------------------------*
* Epilogue                                                             *
*----------------------------------------------------------------------*
      Call qExit('guessorb')
      If(StandAlone) Then
         Call qStat(' ')
         Call FastIO('STATUS')
      End If
      iReturn=0
*----------------------------------------------------------------------*
* Done                                                                 *
*----------------------------------------------------------------------*
      Return
      End
