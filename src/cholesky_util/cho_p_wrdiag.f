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
      SubRoutine Cho_P_WrDiag()
C
C     Purpose: store global diagonal on disk (Parallel only).
C              NB: on exit, initial global diagonal is stored!
C
      use ChoSwp, only: Diag_G
      Implicit None
#include "cholesky.fh"
#include "choglob.fh"
#include "cho_para_info.fh"
#include "stdalloc.fh"
      Real*8, Allocatable:: Diag_L(:)

      If (Cho_Real_Par) Then
         Call mma_allocate(Diag_L,nnBstRT(1),Label='Diag_L')
         Call Cho_IODiag(Diag_L,2)
         Call Cho_P_SyncDiag(Diag_L,1)
         Call Cho_P_IndxSwp()
         Call Cho_IODiag(Diag_G,1)
         Call Cho_P_IndxSwp()
         Call mma_deallocate(Diag_L)
      End If

      End
