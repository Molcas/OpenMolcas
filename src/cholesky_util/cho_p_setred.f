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
      SubRoutine Cho_P_SetRed(Diag,Sync)
C
C     Purpose: set next reduced set after synchronizing the global
C              diagonal (if requested through flag Sync).
C              Global as well as local reduced sets are set.
C              Diag is the local diagonal, whereas ip_Diag_G (stored in
C              choglob.fh) points to the global diagonal in Work.
C              Note that Diag is not referenced if Sync=.False.
C
      Implicit None
      Real*8  Diag(*)
      Logical Sync
#include "cholesky.fh"
#include "WrkSpc.fh"
#include "choglob.fh"
#include "cho_para_info.fh"

      Integer iLoc

      If (Cho_Real_Par) Then

C        Sync global diagonal.
C        ---------------------

         If (Sync) Then
            iLoc = 2
            Call Cho_P_SyncDiag(Diag,iLoc)
         End If

C        Set next global reduced set. The original serial routines
C        are used and so we must trick them by first swapping global
C        and local index arrays (and then swap back, of course).
C        -----------------------------------------------------------

         Call Cho_P_IndxSwp()
         Call Cho_SetRed(Work(ip_Diag_G))
         Call Cho_P_IndxSwp()

C        Set next local reduced set.
C        ---------------------------

         Call Cho_P_SetRed_L()

      Else

         Call Cho_SetRed(Diag)

      End If

      End
