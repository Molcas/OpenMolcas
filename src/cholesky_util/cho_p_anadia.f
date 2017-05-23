************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SubRoutine Cho_P_AnaDia(Diag,Sync,Bin1,Step,NumBin,Full)
C
C     Purpose: analyze global diagonal (histogram). Diag is the local
C              diagonal. If Sync=.True. the global diagonal is
C              synchronized before analysis.
C
      Implicit None
      Real*8  Diag(*)
      Logical Sync
      Real*8  Bin1, Step
      Integer NumBin
      Logical Full
#include "choglob.fh"
#include "cho_para_info.fh"
#include "WrkSpc.fh"

      Integer iLoc

      If (Cho_Real_Par) Then

C        If requested, sync diagonal.
C        ----------------------------

         If (Sync) Then
            iLoc = 2
            Call Cho_P_SyncDiag(Diag,iLoc)
         End If

C        Swap local and global index arrays and use original serial
C        serial to perform analysis of the global diagonal.
C        ----------------------------------------------------------

         Call Cho_P_IndxSwp()
         Call Cho_AnaDia(Work(ip_Diag_G),Bin1,Step,NumBin,Full)
         Call Cho_P_IndxSwp()

      Else

         Call Cho_AnaDia(Diag,Bin1,Step,NumBin,Full)

      End If

      End
