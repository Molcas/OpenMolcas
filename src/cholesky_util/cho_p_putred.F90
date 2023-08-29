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

subroutine Cho_P_PutRed(iPass,iLoc)
!
! Purpose: write global and local reduced set index arrays
!          to disk and set
!          address for next write. iLoc specifies the location in
!          the index arrays to write, and iPass identifies the
!          reduced set (i.e. the integral pass).

use Cholesky, only: Cho_Real_Par, LuRed, LuRed_G, TMISC
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iPass, iLoc
integer(kind=iwp) :: iTmp
real(kind=wp) :: c1, c2, w1, w2

call CWTime(c1,w1)

if (Cho_Real_Par) then

  ! Swap local and global reduced set indices and use original serial
  ! routine to write global index arrays.
  ! -----------------------------------------------------------------

  call Cho_P_IndxSwp()
  iTmp = LuRed
  LuRed = LuRed_G
  call Cho_PutRed(iPass,iLoc)
  LuRed = iTmp
  call Cho_P_IndxSwp()

end if

! Write local index arrays.
! -------------------------

call Cho_PutRed(iPass,iLoc)

call CWTime(c2,w2)
tMisc(1,2) = tMisc(1,2)+c2-c1
tMisc(2,2) = tMisc(2,2)+w2-w1

end subroutine Cho_P_PutRed
