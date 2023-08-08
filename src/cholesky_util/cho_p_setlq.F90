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

subroutine Cho_P_SetLQ()
!
! Purpose: set local qualified indices from known global qualified.
!
! The local indices are stored in the Cholesky module:
!
! nQual_L(iSym)   : #qualified, irrep iSym (=1,2,..,nSym)
! iQuAB_L(iQ,iSym): address of qualified iQ of sym. iSym in current
!                   local reduced set (i.e. reduced set at location
!                   2). Not symmetry reduced, i.e. includes the
!                   offset iiBstR(iSym,2).
! iQL2G(iQ,iSym)  : returns index of the qualified in the global
!                   list.

use Cholesky, only: Cho_Real_Par, iiBstR, iL2G, IndRed, IndRed_G, iQL2G, iQuAB, iQuAB_L, nnBstR, nQual, nQual_L, nSym
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: i, i2, iQ, iQG, iSym, j, k, nQL

if (.not. Cho_Real_Par) return ! not truly parallel...

iQuAB_L(:,:) = 0
iQL2G(:,:) = 0
do iSym=1,nSym
  nQL = 0
  do iQ=1,nQual(iSym)
    iQG = IndRed_G(iQuAB(iQ,iSym),2) ! addr of qual in glob. rs1
    i2 = iiBstR(iSym,2)+nnBstR(iSym,2)
    i = iiBstR(iSym,2)
    do while (i < i2)
      i = i+1
      j = IndRed(i,2) ! addr in local rs1
      k = iL2G(j)     ! addr in global rs1
      if (k == iQG) then ! found qual in local set
        nQL = nQL+1
        iQuAB_L(nQL,iSym) = i
        iQL2G(nQL,iSym) = iQ
        i = i2 ! break while loop
      end if
    end do
  end do
  nQual_L(iSym) = nQL
end do

end subroutine Cho_P_SetLQ
