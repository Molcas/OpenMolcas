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

subroutine Cho_GetLQ(QVec,l_QVec,LstQSP,nQSP)
!
! Purpose: extract elements corresponding to qualified columns from
!          the Cholesky vectors in buffer and/or on disk.

use ChoVecBuf, only: nVec_in_Buf

implicit none
integer l_QVec, nQSP
real*8, target :: QVec(l_Qvec)
integer LstQSP(nQSP)
#include "cholesky.fh"
character*9 SecNam
parameter(SecNam='Cho_GetLQ')
integer iV1(8), nV(8)
integer nTot, iSym
!                                                                      *
!***********************************************************************
!                                                                      *
interface
  subroutine Cho_VecBuf_GetLQ(QVec,l_QVec)
    integer l_QVec
    real*8, target :: QVec(l_QVec)
  end subroutine Cho_VecBuf_GetLQ
end interface
!                                                                      *
!***********************************************************************
!                                                                      *

! Check input.
! ------------

if (nQSP < 1) return
nTot = NumCho(1)
do iSym=2,nSym
  nTot = nTot+NumCho(iSym)
end do
if (nTot < 1) return

nTot = nQual(1)
do iSym=2,nSym
  nTot = nTot+nQual(iSym)
end do
if (nTot < 1) return

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
