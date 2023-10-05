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

subroutine Cho_VecBuf_Copy(Vec,nVec,iSym)
!
! Purpose: copy as many vectors as fit directly into buffer using
!          current reduced set dimension for symmetry iSym. nVec
!          is the number of vectors in array Vec.
!     NB!  It is important that the vector counter NumCho does NOT
!          include the nVec vectors in array Vec.

use Cholesky, only: CHVBUF, ip_CHVBUF_SYM, l_CHVBUF_SYM, nnBstR, NumCho, nVec_in_Buf
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: Vec(*)
integer(kind=iwp), intent(in) :: nVec, iSym
integer(kind=iwp) :: kOff, lCopy, Left, mUsed, nCopy

! Check if there is anything to do at all.
! 1) any buffer allocated?
! 2) any vectors to copy?
! 3) any elements?
! ----------------------------------------

if (l_ChVBuf_Sym(iSym) < 1) return
if (nVec < 1) return
if (nnBstR(iSym,2) < 1) return

! Copy vectors
! 1) if all previous vectors are in the buffer, and
! 2) if there is sufficient free space in the buffer.
! ---------------------------------------------------

if (nVec_in_Buf(iSym) == NumCho(iSym)) then
  mUsed = nnBstR(iSym,2)*nVec_in_Buf(iSym)
  Left = l_ChVBuf_Sym(iSym)-mUsed
  nCopy = min(Left/nnBstR(iSym,2),nVec)
  if (nCopy > 0) then
    lCopy = nnBstR(iSym,2)*nCopy
    kOff = ip_ChVBuf_Sym(iSym)+mUsed-1
    CHVBUF(kOff+1:kOff+lCopy) = Vec(1:lCopy)
    nVec_in_Buf(iSym) = nVec_in_Buf(iSym)+nCopy
  end if
end if

end subroutine Cho_VecBuf_Copy
