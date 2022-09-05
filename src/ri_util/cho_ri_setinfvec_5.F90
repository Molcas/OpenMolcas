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
#if defined (_MOLCAS_MPP_)
      SubRoutine Cho_RI_SetInfVec_5(iVec_Global,iVec_Local,J_s,J_e,iSym)
!
!     Set mapping from local to global vector index (needed in parallel
!     RI gradient code).
!
      use ChoSwp, only: InfVec
      Implicit None
      Integer iVec_Global, iVec_Local, J_s, J_e, iSym
#include "cholesky.fh"

      Integer iOff, nVec, iVec

      iOff = iVec_Global + J_s - 2
      nVec = J_e - J_s + 1
      Do iVec = 1,nVec
         InfVec(iVec_Local-1+iVec,5,iSym) = iOff + iVec
      End Do

      End
#endif
