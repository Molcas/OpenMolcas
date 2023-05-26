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
      SubRoutine Cho_P_SetRed(Diag,Sync)
!
!     Purpose: set next reduced set after synchronizing the global
!              diagonal (if requested through flag Sync).
!              Global as well as local reduced sets are set.
!              Diag is the local diagonal, whereas Diag_G (defined in
!              choswp.f90) points to the global diagonal.
!              Note that Diag is not referenced if Sync=.False.
!
      use ChoSwp, only: Diag_G
      Implicit None
      Real*8  Diag(*)
      Logical Sync
#include "cholesky.fh"
#include "choglob.fh"
#include "cho_para_info.fh"

      Integer iLoc

      If (Cho_Real_Par) Then

!        Sync global diagonal.
!        ---------------------

         If (Sync) Then
            iLoc = 2
            Call Cho_P_SyncDiag(Diag,iLoc)
         End If

!        Set next global reduced set. The original serial routines
!        are used and so we must trick them by first swapping global
!        and local index arrays (and then swap back, of course).
!        -----------------------------------------------------------

         Call Cho_P_IndxSwp()
         Call Cho_SetRed(Diag_G)
         Call Cho_P_IndxSwp()

!        Set next local reduced set.
!        ---------------------------

         Call Cho_P_SetRed_L()

      Else

         Call Cho_SetRed(Diag)

      End If

      End
