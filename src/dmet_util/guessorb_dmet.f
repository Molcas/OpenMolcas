* guessorb $ this file belongs to the Molcas repository $
      subroutine guessorb_dmet(iReturn,StandAlone)
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
*----------------------------------------------------------------------*
*                                                                      *
* Copyright: All rights reserved by Lund University                    *
*                                                                      *
************************************************************************
*      Implicit None
*#include "Molcas.fh"
*#include "commgo.fh"
*----------------------------------------------------------------------*
* Dummy arguments                                                      *
*----------------------------------------------------------------------*
      Integer iReturn
      Logical StandAlone
**----------------------------------------------------------------------*
** Local variables.                                                     *
**----------------------------------------------------------------------*
*      Integer iRC, iUHF
**----------------------------------------------------------------------*
** Prologue                                                             *
**----------------------------------------------------------------------*
      iReturn=0
      Call qEnter('guessorb_dmet')
      Call InitGO(StandAlone)
*      If(StandAlone) Call InpCtl_GuessOrb
**----------------------------------------------------------------------*
** Select method to be used.                                            *
**----------------------------------------------------------------------*
**VB*   Call cre_gsswfn
*      Call FckByInt(iRC,StandAlone)
**     If(iRC.eq.0) GoTo 999
**VB*  If(.true.) GoTo 999
**VB*  If(nSym.eq.1) Then
**VB*     Call Fmod1n(StandAlone)
**VB*  Else
**VB*     Call Fmod1s(StandAlone)
**VB*  End If
*999   Continue
*      Call cls_gsswfn
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
