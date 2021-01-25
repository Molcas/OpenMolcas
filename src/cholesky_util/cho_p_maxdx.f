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
      SubRoutine Cho_P_MaxDX(Diag,Sync,Dmax)
C
C     Purpose: get max. diagonal elements in each sym. block,
C              qualified diagonals excluded.
C
      use ChoSwp, only: Diag_G
      Implicit None
      Real*8  Diag(*)
      Logical Sync
      Real*8  Dmax(*)
#include "cho_para_info.fh"
#include "choglob.fh"

      Integer iLoc

      If (Cho_Real_Par) Then
         If (Sync) Then
            iLoc = 2
            Call Cho_P_SyncDiag(Diag,iLoc)
         End If
         Call Cho_P_IndxSwp()
         Call Cho_MaxDX(Diag_G,Dmax)
         Call Cho_P_IndxSwp()
      Else
         Call Cho_MaxDX(Diag,Dmax)
      End If

      End
