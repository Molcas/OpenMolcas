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
      SubRoutine Cho_P_IndxSwp()
C
C     Purpose: swap global and local reduced set index arrays.
C              Note that arrays InfRed and InfVec are swapped as well.
C
C     NB: this procedure is inexpensive, as we are merely swapping
C         pointers, not actual data (except for the statically allocated
C         index arrays which amount to swapping 51 integers in total).
C
      use ChoSwp, only: nnBstRSh, nnBstRSh_G, pTemp3
      use ChoSwp, only: iiBstRSh, iiBstRSh_G
      use ChoSwp, only: IndRSh, IndRSh_G, pTemp1
      use ChoSwp, only: InfRed, InfRed_G
      Implicit None
#include "cholesky.fh"
#include "choptr.fh"
#include "choglob.fh"
#include "WrkSpc.fh"

      Integer iTmp, N

      iTmp = nnShl_G
      nnShl_G = nnShl
      nnShl = iTmp

      iTmp = mmBstRT_G
      mmBstRT_G = mmBstRT
      mmBstRT = iTmp

      N = 8*3
      Call iSwap(N,iiBstR_G,1,iiBstR,1)
      Call iSwap(N,nnBstR_G,1,nnBstR,1)
      Call iSwap(3,nnBstRT_G,1,nnBstRT,1)

      pTemp1 => InfRed_G
      InfRed_G => InfRed
      InfRed => pTemp1

      iTmp = ip_InfVec_G
      ip_InfVec_G = ip_InfVec
      ip_InfVec = iTmp
      iTmp = l_InfVec_G
      l_InfVec_G = l_InfVec
      l_InfVec = iTmp

      pTemp3 => iiBstRSh_G
      iiBstRSh_G => iiBstRSh
      iiBstRSh => pTemp3

      pTemp3 => nnBstRSh_G
      nnBstRSh_G => nnBstRSh
      nnBstRSh => pTemp3

      iTmp = ip_IndRed_G
      ip_IndRed_G = ip_IndRed
      ip_IndRed = iTmp
      iTmp = l_IndRed_G
      l_IndRed_G = l_IndRed
      l_IndRed = iTmp

      pTemp1 => IndRSh_G
      IndRSh_G => IndRSh
      IndRSh => pTemp1

      End
