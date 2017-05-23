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
      SubRoutine Cho_P_PrtDia(Diag,Sync,iSyLst,nSyLst,iLoc)
C
C     Purpose: print global diagonal. Diag is the local diagonal and
C              if Sync=.True. the global diagonal is
C              synchronized before printing. Array iSyLst(nSyLst)
C              specifies which symmetry blocks to print, and iLoc points
C              to the memory location of the reduced set index arrays to
C              use for printing (and synchronizing, if requested).
C
      Implicit None
      Real*8  Diag(*)
      Logical Sync
      Integer nSyLst
      Integer iSyLst(nSyLst)
      Integer iLoc
#include "cho_para_info.fh"
#include "choglob.fh"
#include "WrkSpc.fh"

      If (Cho_Real_Par) Then

C        Sync diagonal if requested.
C        ---------------------------

         If (Sync) Then
            Call Cho_P_SyncDiag(Diag,iLoc)
         End If

C        Swap local and global index arrays and use original serial routine
C        to print diagonal.
C        ------------------------------------------------------------------

         Call Cho_P_IndxSwp()
         Call Cho_PrtDia(Work(ip_Diag_G),iSyLst,nSyLst,iLoc)
         Call Cho_P_IndxSwp()

      Else

         Call Cho_PrtDia(Diag,iSyLst,nSyLst,iLoc)

      End If

      End
