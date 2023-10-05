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
! Copyright (C) 2012, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_P_UpdateBookmarks(iRS)
!
! Thomas Bondo Pedersen, August 2012.
!
! Update bookmarks for reduced set (integral pass) iRS:
!    - integral accuracy (max diagonal)
!    - number of Cholesky vectors
!
! Note: it is assumed that array DiaMax and number of Cholesky
! vectors are properly updated before calling this routine.

use Cholesky, only: BkmThr, BkmVec, Cho_Real_Par, DiaMaxT, nCol_BkmThr, nCol_BkmVec, nSym, NumCho, NumCho_G
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: iRS

if ((.not. allocated(BkmVec)) .or. (.not. allocated(BkmThr))) return

if (Cho_Real_Par) then
  BkmVec(1:nSym,iRS) = NumCho_G(1:nSym)
else
  BkmVec(1:nSym,iRS) = NumCho(1:nSym)
end if
BkmThr(1:nSym,iRS) = DiaMaxT(1:nSym)
nCol_BkmVec = nCol_BkmVec+1
nCol_BkmThr = nCol_BkmThr+1

end subroutine Cho_P_UpdateBookmarks
