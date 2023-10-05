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

subroutine Cho_P_SyncDiag(Diag,iLoc)
!
! Purpose: synchronize global diagonal. Diag is the local diagonal
!          and iLoc tells which memory location to use for reduced
!          set index arrays.

use Cholesky, only: Cho_Real_Par, Diag_G, iL2G, IndRed, nnBstRT, nnBstRT_G, TMISC
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: Diag(*)
integer(kind=iwp), intent(in) :: iLoc
integer(kind=iwp) :: i, j
real(kind=wp) :: c1, c2, w1, w2

! Skip if serial run.
! -------------------

if (.not. Cho_Real_Par) return
call CWTime(c1,w1)

! Zero all entries in global diagonal.
! ------------------------------------

Diag_G(1:nnBstRT_G(1)) = Zero

! Copy elements from local to global diagonal.
! --------------------------------------------

if (iLoc == 1) then
  do i=1,nnBstRT(1)
    Diag_G(iL2G(i)) = Diag(i)
  end do
else
  do j=1,nnBstRT(iLoc)
    i = IndRed(j,iLoc)
    Diag_G(iL2G(i)) = Diag(i)
  end do
end if

! Synchronize global diagonal.
! ----------------------------

call Cho_GADGop(Diag_G,nnBstRT_G(1),'+')

call CWTime(c2,w2)
tMisc(1,4) = tMisc(1,4)+c2-c1
tMisc(2,4) = tMisc(2,4)+w2-w1

end subroutine Cho_P_SyncDiag
