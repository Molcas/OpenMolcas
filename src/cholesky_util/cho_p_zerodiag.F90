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

subroutine Cho_P_ZeroDiag(Diag,iSym,iABG)
!
! Purpose: zero diagonal element iABG (in global diagonal, rs1).
!          For serial runs, this is trivial. For parallel runs, we
!          need first to figure out if the treated diagonal element
!          is in fact present among the qualified in the local
!          diagonal.
!
! NB! If you wish to test the entire local diagonal (i.e. not just
!     the qualified), use Cho_P_ZeroDiag_Rst instead.

use Cholesky, only: Cho_Real_Par, iL2G, IndRed, iQuAB_L, nQual_L
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: Diag(*)
integer(kind=iwp), intent(in) :: iSym, iABG
integer(kind=iwp) :: iAB, iQ, jAB, kAB

if (Cho_Real_Par) then
  do iQ=1,nQual_L(iSym)
    iAB = iQuAB_L(iQ,iSym) ! addr in local current rs
    jAB = IndRed(iAB,2)    ! addr in local rs1
    kAB = iL2G(jAB)        ! addr in global rs1
    if (kAB == iABG) then  ! found...
      Diag(jAB) = Zero     ! now zero local diagonal elm.
      return
    end if
  end do
else
  Diag(iABG) = Zero
end if

end subroutine Cho_P_ZeroDiag
