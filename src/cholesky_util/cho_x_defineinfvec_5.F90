!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************

subroutine Cho_X_DefineInfVec_5(isDF)
!
! Purpose: Trivial definition of location 5 of InfVec:
!          InfVec(i,5,iSym) = i
!          The routine does nothing in case of parallel DF.

use Para_Info, only: Is_Real_Par
use Cholesky, only: InfVec, nSym, NumCho
use Definitions, only: iwp

implicit none
logical(kind=iwp), intent(in) :: isDF
integer(kind=iwp) :: i, iSym
logical(kind=iwp) :: doDefine

! Define in case of
! 1) serial Cholesky
! 2) serial DF
! 3) parallel Cholesky
! Do NOT define for parallel DF.
doDefine = .not. Is_Real_Par() .or. (Is_Real_Par() .and. (.not. isDF))
if (doDefine) then
  do iSym=1,nSym
    do i=1,NumCho(iSym)
      InfVec(i,5,iSym) = i
    end do
  end do
end if

end subroutine Cho_X_DefineInfVec_5
