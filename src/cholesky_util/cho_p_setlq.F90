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
! The local indices are stored in ChoArr.f90 and ChoSwp.f90:
!
! nQual_L(iSym)   : #qualified, irrep iSym (=1,2,..,nSym)
! iQuAB_L(iQ,iSym): address of qualified iQ of sym. iSym in current
!                   local reduced set (i.e. reduced set at location
!                   2). Not symmetry reduced, i.e. includes the
!                   offset iiBstR(iSym,2).
! iQL2G(iQ,iSym)  : returns index of the qualified in the global
!                   list.

use ChoSwp, only: iQuAB, iQuAB_L, IndRed, IndRed_G
use ChoArr, only: iL2G, iQL2G, nQual_L

implicit none
#include "cholesky.fh"
#include "choglob.fh"
#include "cho_para_info.fh"
integer iSym, nQL, iQ, iQG, i2, i, j, k

if (.not. Cho_Real_Par) return ! not truly parallel...

call iZero(iQuAB_L,size(iQuAB_L))
call iZero(iQL2G,size(iQL2G))
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
