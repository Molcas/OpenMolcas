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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine Drv2El_LocalDF()
C
C     Thomas Bondo Pedersen, June 2010.
C
C     Purpose: Compute fitting coefficients for Local Density Fitting.
C
      Implicit None
#include "localdf.fh"
#include "localdf_print.fh"

      Character*14 SecNam
      Parameter (SecNam='Drv2El_LocalDF')

      Logical DoPairs
      Logical Verbose
      Logical TermInts

      Integer irc, i

      Real*8 tCPU0, tCPU1, tC0, tC1
      Real*8 tWall0, tWall1, tW0, tW1

      Integer  iPrintLevel
      External iPrintLevel

      Call qEnter('2El_LDF')
*

      Call StatusLine('Seward: ','local density fitting')
#if defined (_DEBUG_)
      Write(6,'(A,A)') '>>> Enter ',SecNam
#endif
      Call ThisIsRestrictedCode('Thomas Bondo Pedersen',
     &                          'Local Density Fitting',.false.)

C=============
C     Preamble
C=============

      ! Set print level
      iPrint=iPrintLevel(-1)
      ! Start overall timing
      If (iPrint.ge.Inf_OverallTiming) Then
         Call CWTime(tCPU0,tWall0)
      End If
      ! Init return code
      irc=0
      ! Set run mode to "internal"
      LDF_Run_Mode=LDF_RUN_INTERNAL
      ! Set LuPri in Cholesky (enables use of printing routines in
      ! Cholesky utility)
      Call LDF_Cho_SetLuPri(6)
      ! Print header
      If (iPrint.ge.Inf_Min) Then
         Write(6,'(//,80A1)') ('*',i=1,80)
         Write(6,'(A1,78X,A1)') ('*',i=1,2)
         Write(6,'(A1,10X,A,10X,A1)')
     &   '*',
     &   'Local Density Fitting: Calculation of Fitting Coefficients',
     &   '*'
         Write(6,'(A1,78X,A1)') ('*',i=1,2)
         Write(6,'(80A1)') ('*',i=1,80)
         If (LDF2) Then
            Write(6,'(A)')
     &   'Inclusion of two-center auxiliary functions...             ON'
         Else
            Write(6,'(A)')
     &   'Inclusion of two-center auxiliary functions...            OFF'
         End If
         If (LDF_Constraint.eq.-1) Then
            Write(6,'(A)')
     &   'Constraint....................................           None'
         Else If (LDF_Constraint.eq.0) Then
            Write(6,'(A)')
     &   'Constraint....................................         Charge'
         Else
            Write(6,'(A)')
     &   'Constraint....................................        Unknown'
         End If
         Write(6,'(A,1P,D15.6)')
     &   'Target Accuracy...............................',Thr_Accuracy
         Write(6,'(A,1P,D15.6)')
     &   'Prescreening Threshold........................',Thr_Prescreen
         Write(6,'(A,5X,I10)')
     &   'Print Level...................................',iPrint
         Call xFlush(6)
      End If

C====================
C     Initialization.
C====================

      If (iPrint.ge.Inf_DetailedTiming) Then
         Write(6,'(/,A)')
     &   '***** Starting LDF initialization *****'
         Call xFlush(6)
         Call CWTime(tC0,tW0)
      End If
      DoPairs=.True.
      Verbose=iPrint.ge.Inf_Init
      Call LDF_Init(DoPairs,Verbose,irc)
      If (irc.ne.0) Then
         Write(6,'(A,A,I8)')
     &   SecNam,': LDF_Init returned code',irc
         Call LDF_Quit(irc)
      End If
      If (iPrint.ge.Inf_DetailedTiming) Then
         Call CWTime(tC1,tW1)
         Call Cho_PrtTim('LDF initialization',tC1,tC0,tW1,tW0,1)
      End If

C==================================
C     Compute fitting coefficients.
C==================================

      If (iPrint.ge.Inf_DetailedTiming) Then
         Write(6,'(/,A)')
     &   '***** Starting LDF fitting *****'
         Call xFlush(6)
         Call CWTime(tC0,tW0)
      End If
      Call LDF_ComputeFittingCoefficients(irc)
      If (irc.ne.0) Then
         Write(6,'(A,A,I8)')
     &   SecNam,': LDF_ComputeFittingCoefficients returned code',irc
         Call LDF_Quit(irc)
      End If
      If (iPrint.ge.Inf_DetailedTiming) Then
         Call CWTime(tC1,tW1)
         Call Cho_PrtTim('LDF fitting',tC1,tC0,tW1,tW0,1)
      End If

C==================
C     Finalization.
C==================

      If (iPrint.ge.Inf_DetailedTiming) Then
         Write(6,'(/,A)')
     &   '***** Starting LDF finalization *****'
         Call xFlush(6)
         Call CWTime(tC0,tW0)
      End If
      Call Put_dScalar('Cholesky Threshold',Thr_Accuracy)
      TermInts=.True.
      Call LDF_Final(TermInts,irc)
      If (irc.ne.0) Then
         Write(6,'(A,A,I8)')
     &   SecNam,': LDF_Final returned code',irc
         Call LDF_Quit(irc)
      End If
      If (iPrint.ge.Inf_DetailedTiming) Then
         Call CWTime(tC1,tW1)
         Call Cho_PrtTim('LDF finalization',tC1,tC0,tW1,tW0,1)
      End If
      Call Free_iSD()

      ! Print overall timing
      If (iPrint.ge.Inf_OverallTiming) Then
         Call CWTime(tCPU1,tWall1)
         Call Cho_PrtTim('Local Density Fitting',
     &                   tCPU1,tCPU0,tWall1,tWall0,1)
      End If

#if defined (_DEBUG_)
      Write(6,'(A,A)') '>>> Exit ',SecNam
#endif
      Call qExit('2El_LDF')
      Return
      End
