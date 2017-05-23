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
* Copyright (C) 2008, Francesco Aquilante                              *
************************************************************************
      SubRoutine ChoMP2_FNO(irc,D_ab,D_ii,EOcc,EVir,Sorted,DelOrig)
C
C     F. Aquilante, Geneva May 2008  (snick in Pedersen's code)
C
C
      Implicit None
      Integer irc
      Real*8  D_ab(*), D_ii(*)
      Real*8  EOcc(*), EVir(*)
      Logical Sorted, DelOrig
#include "chomp2.fh"
#include "WrkSpc.fh"

      Character*3  ThisNm
      Character*10 SecNam
      Parameter (SecNam = 'ChoMP2_FNO', ThisNm = 'FNO')

      Integer ipWrk, lWrk

      Call qEnter(ThisNm)
      irc = 0

      Call GetMem('GetMax','Max ','Real',ipWrk,lWrk)
      Call GetMem('GetMax','Allo','Real',ipWrk,lWrk)

      If (Sorted) Then
         Call ChoMP2_fno_Srt(irc,DelOrig,D_ab,D_ii,EOcc,EVir,
     &                           Work(ipWrk),lWrk)
         If (irc .ne. 0) Then
            Write(6,*) SecNam,': ChoMP2_fno_Srt returned ',irc
            Go To 1 ! exit
         End If
      Else
         If (nBatch .eq. 1) Then
            Call ChoMP2_fno_Fll(irc,DelOrig,D_ab,D_ii,EOcc,EVir,
     &                              Work(ipWrk),lWrk)
            If (irc .ne. 0) Then
               Write(6,*) SecNam,': ChoMP2_fno_Fll returned ',irc
               Go To 1 ! exit
            End If
         Else
            Call ChoMP2_fno_Org(irc,DelOrig,D_ab,D_ii,EOcc,EVir,
     &                              Work(ipWrk),lWrk)
            If (irc .ne. 0) Then
               Write(6,*) SecNam,': ChoMP2_fno_Org returned ',irc
               Go To 1 ! exit
            End If
         End If
      End If

    1 Call GetMem('GetMax','Free','Real',ipWrk,lWrk)
      Call qExit(ThisNm)
      End
