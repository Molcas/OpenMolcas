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
* Copyright (C) 2004, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine ChoMP2_Energy(irc,EMP2,EOcc,EVir,Sorted,DelOrig)
C
C     Thomas Bondo Pedersen, Dec. 2004.
C
C     Purpose: compute MP2 energy correction from MO Cholesky vectors,
C              constructing (ai|bj) integrals on the fly. Flag Sorted
C              refers to whether or not the MO vectors have been sorted
C              into the sizes of the batches over occupied orbitals.
C
      Implicit None
      Real*8  EMP2
      Real*8  EOcc(*), EVir(*)
      Integer irc
      Logical Sorted, DelOrig
#include "chomp2.fh"
#include "chomp2_cfg.fh"
#include "WrkSpc.fh"

      Character*6  ThisNm
      Character*13 SecNam
      Parameter (SecNam = 'ChoMP2_Energy', ThisNm = 'Energy')

      Integer ipWrk, lWrk

      irc = 0

      Call GetMem('GetMax','Max ','Real',ipWrk,lWrk)
      Call GetMem('GetMax','Allo','Real',ipWrk,lWrk)

      If (Sorted) Then
         Call ChoMP2_Energy_Srt(irc,DelOrig,EMP2,EOcc,EVir,
     &                          Work(ipWrk),lWrk)
         If (irc .ne. 0) Then
            Write(6,*) SecNam,': ChoMP2_Energy_Srt returned ',irc
            Go To 1 ! exit
         End If
      Else
         If (nBatch .eq. 1) Then
            Call ChoMP2_Energy_Fll(irc,DelOrig,EMP2,EOcc,EVir,
     &                             Work(ipWrk),lWrk)
            If (irc .ne. 0) Then
               Write(6,*) SecNam,': ChoMP2_Energy_Fll returned ',irc
               Go To 1 ! exit
            End If
         Else
            Call ChoMP2_Energy_Org(irc,DelOrig,EMP2,EOcc,EVir,
     &                             Work(ipWrk),lWrk)
            If (irc .ne. 0) Then
               Write(6,*) SecNam,': ChoMP2_Energy_Org returned ',irc
               Go To 1 ! exit
            End If
         End If
      End If

    1 Call GetMem('GetMax','Free','Real',ipWrk,lWrk)
      End
