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

subroutine Cho_P_SetVecInf(nVec,iSym,iPass)
!
! Purpose: set global and local info for vectors.

use Cholesky, only: Cho_Real_Par, iQuAB, IndRed, NumCho, NumCho_G
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nVec, iSym, iPass
integer(kind=iwp) :: iV, iVec, iAB
integer(kind=iwp), external :: Cho_P_IndxParentDiag

if (Cho_Real_Par) then
  ! Set global vector information (by swapping index arrays)
  call Cho_P_IndxSwp()
  do iV=1,nVec
    iVec = NumCho_G(iSym)+iV
    iAB = IndRed(iQuAB(iV,iSym),2)
    call Cho_SetVecInf(iVec,iSym,iAB,iPass,2)
  end do
  call Cho_P_IndxSwp()
  ! Set local vector information
  do iV=1,nVec
    iVec = NumCho_G(iSym)+iV
    iAB = Cho_P_IndxParentDiag(iV,iSym)
    call Cho_SetVecInf(iVec,iSym,iAB,iPass,2)
  end do
else
  do iV=1,nVec
    iVec = NumCho(iSym)+iV
    iAB = IndRed(iQuAB(iV,iSym),2)
    call Cho_SetVecInf(iVec,iSym,iAB,iPass,2)
  end do
end if

end subroutine Cho_P_SetVecInf
