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
      SubRoutine Cho_P_SetPass(Diag,Sync,DiaSh,iSySh,iLoc,Conv,nPotSh)
C
C     Purpose: check convergence and, if not converged, set up next
C              next integral pass
C              Diag is the local diagonal; the global diagonal is
C              synchronized if Sync=.True. Note that DiaSh and iSySh
C              must be allocated with dimension nnShl_G (global number
C              of shell pairs in 1st reduced set). iLoc is the location
C              to use in the reduced set index arrays. On exit,
C              Conv=.True. if converged and nPotSh is the number of
C              shell pairs whose max. diagonal element is larger than
C              the decomposition threshold.
C
      use ChoSwp, only: Diag_G
      Implicit None
      Real*8  Diag(*)
      Logical Sync, Conv
      Real*8  DiaSh(*)
      Integer iSySh(*)
      Integer iLoc, nPotSh
#include "cho_para_info.fh"
#include "choglob.fh"

      If (Cho_Real_Par) Then

C        Sync diagonal if requested.
C        ---------------------------

         If (Sync) Then
            Call Cho_P_SyncDiag(Diag,iLoc)
         End If

C        Swap local and global index arrays and set next integral pass
C        original serial routine.
C        -------------------------------------------------------------

         Call Cho_P_IndxSwp()
         Call Cho_SetPass(Diag_G,DiaSh,iSySh,iLoc,Conv,nPotSh)
         Call Cho_P_IndxSwp()

      Else

         Call Cho_SetPass(Diag,DiaSh,iSySh,iLoc,Conv,nPotSh)

      End If

      End
