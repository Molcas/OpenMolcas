!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      SubRoutine Cho_P_Qualify(Diag,Sync,iShlAB,iSyMax,Mem,Full)
!
!     Purpose: qualify global diagonal elements for decomposition in
!              current reduced set.
!              iShlAB is the shell pair from which to qualify.
!              iSyMax is the symmetry block to which the largest
!              diagonal belongs.
!              Mem is the (total!) max. allowed
!              memory for storing the qualified columns.
!              If no more columns can be qualified on exit,
!              FULL=.true. is returned. On entry, Diag is the local
!              diagonal. The global diagonal is synchronized if
!              Sync=.True. on entry.
!
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

!        Sync global diagonal (if requested).
!        ------------------------------------

         If (Sync) Then
            iLoc = 2
            Call Cho_P_SyncDiag(Diag,iLoc)
         End If

!        Swap local and global index arrays and use the original serial
!        routines for qualifying diagonals.
!        --------------------------------------------------------------
!-TODO/FIXME: the memory constraint is then based on the global
!             dimension, thus ensuring that all nodes find the same
!             qualified (and this is absolutely essential!). However,
!             one might find a more clever way of ensuring this....

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
