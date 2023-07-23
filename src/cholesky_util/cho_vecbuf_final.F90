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

subroutine Cho_VecBuf_Final()
!
! Purpose: deallocate and finalize vector buffer.

use ChoVecBuf, only: CHVBFI, CHVBUF, ip_CHVBFI_SYM, ip_CHVBUF_SYM, l_CHVBFI_SYM, l_CHVBUF_SYM, nVec_in_Buf
use stdalloc, only: mma_deallocate

implicit none
#include "cholesky.fh"

if (allocated(CHVBUF)) call mma_deallocate(CHVBUF)
if (allocated(CHVBFI)) call mma_deallocate(CHVBFI)

call iZero(ip_ChVBuf_Sym,nSym)
call iZero(l_ChVBuf_Sym,nSym)
call iZero(ip_ChVBFI_Sym,nSym)
call iZero(l_ChVBFI_Sym,nSym)
call iZero(nVec_in_Buf,nSym)

end subroutine Cho_VecBuf_Final
