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

use Cholesky, only: BkmThr, BkmVec, Cho_Real_Par, DiaMaxT, MaxRed, nCol_BkmThr, nCol_BkmVec, nSym, NumCho, NumCho_G
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iRS

if ((.not. allocated(BkmVec)) .or. (.not. allocated(BkmThr))) return

if (Cho_Real_Par) then
  call Cho_UpdateBookmarks(iRS,nSym,MaxRed,NumCho_G,DiaMaxT,BkmVec,BkmThr)
else
  call Cho_UpdateBookmarks(iRS,nSym,MaxRed,NumCho,DiaMaxT,BkmVec,BkmThr)
end if
nCol_BkmVec = nCol_BkmVec+1
nCol_BkmThr = nCol_BkmThr+1

end subroutine Cho_P_UpdateBookmarks
