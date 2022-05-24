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
      SubRoutine Cho_P_Qualify(Diag,Sync,iShlAB,iSyMax,Mem,Full)
C
C     Purpose: qualify global diagonal elements for decomposition in
C              current reduced set.
C              iShlAB is the shell pair from which to qualify.
C              iSyMax is the symmetry block to which the largest
C              diagonal belongs.
C              Mem is the (total!) max. allowed
C              memory for storing the qualified columns.
C              If no more columns can be qualified on exit,
C              FULL=.true. is returned. On entry, Diag is the local
C              diagonal. The global diagonal is synchronized if
C              Sync=.True. on entry.
C
      use ChoSwp, only: Diag_G
      Implicit None
      Real*8  Diag(*)
      Logical Sync, Full
      Integer iShlAB, iSyMax, Mem
#include "cholesky.fh"
#include "cho_para_info.fh"
#include "choglob.fh"

      Integer iLoc

      Real*8 c1, c2, w1, w2

      Call Cho_Timer(c1,w1)

      If (Cho_Real_Par) Then

C        Sync global diagonal (if requested).
C        ------------------------------------

         If (Sync) Then
            iLoc = 2
            Call Cho_P_SyncDiag(Diag,iLoc)
         End If

C        Swap local and global index arrays and use the original serial
C        routines for qualifying diagonals.
C        --------------------------------------------------------------
C-TODO/FIXME: the memory constraint is then based on the global
C             dimension, thus ensuring that all nodes find the same
C             qualified (and this is absolutely essential!). However,
C             one might find a more clever way of ensuring this....

         Call Cho_P_IndxSwp()
         Call Cho_Qualify(Diag_G,iShlAB,iSyMax,Mem,Full)
         Call Cho_P_IndxSwp()

      Else

         Call Cho_Qualify(Diag,iShlAB,iSyMax,Mem,Full)

      End If

      Call Cho_Timer(c2,w2)
      tMisc(1,1)=tMisc(1,1)+c2-c1
      tMisc(2,1)=tMisc(2,1)+w2-w1

      End
