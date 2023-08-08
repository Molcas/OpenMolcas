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

use Cholesky, only: CHVBFI, CHVBUF, ip_CHVBFI_SYM, ip_CHVBUF_SYM, l_CHVBFI_SYM, l_CHVBUF_SYM, nSym, nVec_in_Buf
use stdalloc, only: mma_deallocate

implicit none

if (allocated(CHVBUF)) call mma_deallocate(CHVBUF)
if (allocated(CHVBFI)) call mma_deallocate(CHVBFI)

ip_ChVBuf_Sym(1:nSym) = 0
l_ChVBuf_Sym(1:nSym) = 0
ip_ChVBFI_Sym(1:nSym) = 0
l_ChVBFI_Sym(1:nSym) = 0
nVec_in_Buf(1:nSym) = 0

end subroutine Cho_VecBuf_Final
