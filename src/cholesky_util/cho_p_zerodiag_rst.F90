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

subroutine Cho_P_ZeroDiag_Rst(Diag,iSym,iABG)
!
! Purpose: zero diagonal element iABG (in global diagonal, rs1).
!          For serial runs, this is trivial. For parallel runs, we
!          need first to figure out if the treated diagonal element
!          is in fact present in the local diagonal.

use ChoSwp, only: IndRed
use ChoArr, only: iL2G

implicit none
real*8 Diag(*)
integer iSym, iABG
#include "cho_para_info.fh"
#include "cholesky.fh"
#include "choglob.fh"
integer iAB1, iAB2, iAB, jAB, kAB

if (Cho_Real_Par) then
  iAB1 = iiBstR(iSym,2)+1
  iAB2 = iAB1+nnBstR(iSym,2)-1
  do iAB=iAB1,iAB2
    jAB = IndRed(iAB,2)   ! addr in local rs1
    kAB = iL2G(jAB)       ! addr in global rs1
    if (kAB == iABG) then ! found...
      Diag(jAB) = 0.0d0   ! now zero local diagonal elm.
      return
    end if
  end do
else
  Diag(iABG) = 0.0d0
end if

end subroutine Cho_P_ZeroDiag_Rst
