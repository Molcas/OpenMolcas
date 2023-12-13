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

subroutine Cho_GetLQ(QVec,l_QVec,LstQSP,nQSP)
!
! Purpose: extract elements corresponding to qualified columns from
!          the Cholesky vectors in buffer and/or on disk.

use Cholesky, only: nQual, nSym, NumCho, nVec_in_Buf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: l_QVec, nQSP, LstQSP(nQSP)
real(kind=wp), target, intent(out) :: QVec(l_Qvec)
integer(kind=iwp) :: iSym, iV1(8), nV(8)

! Check input.
! ------------

if (nQSP < 1) return
if (sum(NumCho(1:nSym)) < 1) return

if (sum(nQual(1:nSym)) < 1) return

! Extract from vectors in buffer.
! -------------------------------

call Cho_VecBuf_GetLQ(QVec,l_QVec)

! Extract from vectors on disk.
! -----------------------------

do iSym=1,nSym
  iV1(iSym) = nVec_in_Buf(iSym)+1
  nV(iSym) = NumCho(iSym)-nVec_in_Buf(iSym)
end do
call Cho_VecDsk_GetLQ(QVec,l_QVec,LstQSP,nQSP,iV1,nV,nSym)

end subroutine Cho_GetLQ
