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

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

subroutine Cho_VecBuf_GetLQ(QVec,l_QVec)
!
! Purpose: extract elements corresponding to qualified diagonals
!          from vectors in buffer.

use Cholesky, only: CHVBUF, iiBstR, ip_CHVBUF_SYM, iQuAB, nnBstR, nQual, nSym, nVec_in_Buf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: l_QVec
real(kind=wp), target, intent(out) :: QVec(l_QVec)
integer(kind=iwp) :: iAB, iE, iQ, iS, iSym, iVec, kOffQ, lCol, lRow, nVecTot(8)
real(kind=wp), pointer :: BVec(:,:), Q(:,:)

! Check if there is any buffer at all.
! ------------------------------------

if (.not. allocated(CHVBUF)) return

! Extract from each symmetry block.
! ---------------------------------

call Cho_P_GetGV(nVecTot,nSym)

kOffQ = 0
do iSym=1,nSym

  lRow = nnBstR(iSym,2)
  lCol = nVec_in_Buf(iSym)
  iS = ip_ChVBuf_Sym(iSym)
  iE = iS-1+lRow*lCol
  BVec(1:lRow,1:lCol) => CHVBUF(iS:iE)

  lRow = nQual(iSym)
  lCol = nVec_in_Buf(iSym)
  iS = kOffQ+1
  iE = iS-1+lRow*lCol
  Q(1:lRow,1:lCol) => QVec(iS:iE)

  if (nQual(iSym) > 0) then
    do iVec=1,nVec_in_Buf(iSym)
      do iQ=1,nQual(iSym)
        iAB = iQuAB(iQ,iSym)-iiBstR(iSym,2)
        Q(iQ,iVec) = BVec(iAB,iVec)
      end do
    end do
    kOffQ = kOffQ+nQual(iSym)*nVecTot(iSym)
  end if
end do

nullify(BVec)
nullify(Q)

end subroutine Cho_VecBuf_GetLQ
