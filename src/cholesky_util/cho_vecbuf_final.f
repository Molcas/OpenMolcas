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
      SubRoutine Cho_VecBuf_Final()
!
!     Purpose: deallocate and finalize vector buffer.
!
      use ChoVecBuf, only: CHVBUF, ip_CHVBUF_SYM, l_CHVBUF_SYM,         &
     &                     CHVBFI, ip_CHVBFI_SYM, l_CHVBFI_SYM,         &
     &                     nVec_in_Buf
      use stdalloc
      Implicit None
#include "cholesky.fh"

      If (Allocated(CHVBUF)) Call mma_deallocate(CHVBUF)
      If (Allocated(CHVBFI)) Call mma_deallocate(CHVBFI)

      Call iZero(ip_ChVBuf_Sym,nSym)
      Call iZero(l_ChVBuf_Sym,nSym)
      Call iZero(ip_ChVBFI_Sym,nSym)
      Call iZero(l_ChVBFI_Sym,nSym)
      Call iZero(nVec_in_Buf,nSym)

      End
