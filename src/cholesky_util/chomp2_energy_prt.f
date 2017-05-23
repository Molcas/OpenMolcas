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
* Copyright (C) 2005, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine ChoMP2_Energy_Prt(Caller,Job,iBatch)
C
C     Thomas Bondo Pedersen, March 2005.
C
C     Purpose: print progress reports for the MP2 energy evaluation.
C
      Implicit None
      Character*17 Caller
      Integer      Job, iBatch
#include "chomp2.fh"

      Real*8 CME_Time
      Common / CMETim / CME_Time(2,2)

      Character*17 SecNam
      Parameter (SecNam = 'ChoMP2_Energy_Prt')

      Character*10 ThisNam
      Parameter (ThisNam = 'Energy_Prt')

      Real*8 CPU, Wall, Ratio

      If (Job .eq. 0) Then

         Call Cho_dZero(CME_Time,2*2)

         Write(6,'(/,4X,A,/,4X,A)')
     &   'Evaluation of MP2 energy correction',
     &   '==================================='
         Write(6,'(4X,A,A)')
     &   'Evaluator: ',Caller

         Write(6,'(/,4X,A,/,4X,A,/,4X,A)')
     &   'Batch      CPU       Wall    Ratio',
     &   ' No.     seconds    seconds',
     &   '----------------------------------'

         Call xFlush(6)

      Else If (Job .eq. 1) Then

         Call CWTime(CME_Time(1,1),CME_Time(2,1))

         Call xFlush(6)

      Else If (Job .eq. 2) Then

         Call CWTime(CME_Time(1,2),CME_Time(2,2))
         CPU   = CME_Time(1,2) - CME_Time(1,1)
         Wall  = CME_Time(2,2) - CME_Time(2,1)
         If (Abs(Wall) .lt. 1.0D-8) Then
            If (Abs(CPU) .lt. 1.0D-8) Then
               Ratio = 1.0D0
            Else
               Ratio = 1.0D15
            End If
         Else
            Ratio = CPU/Wall
         End If
         Write(6,'(I9,2(1X,F10.2),1X,F6.3)') iBatch,CPU,Wall,Ratio

         Call xFlush(6)

      Else If (Job .eq. 3) Then

         Write(6,'(4X,A)')
     &   '----------------------------------'

         Call xFlush(6)

      Else

         Call qEnter(ThisNam)
         Call ChoMP2_Quit(SecNam,
     &                    'Input parameter "Job" is out of range',' ')

      End If

      End
