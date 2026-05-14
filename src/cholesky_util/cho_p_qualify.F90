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

subroutine Cho_P_Qualify(Diag,Sync,iShlAB,iSyMax,Mem,Full)
!
! Purpose: qualify global diagonal elements for decomposition in
!          current reduced set.
!          iShlAB is the shell pair from which to qualify.
!          iSyMax is the symmetry block to which the largest
!          diagonal belongs.
!          Mem is the (total!) max. allowed
!          memory for storing the qualified columns.
!          If no more columns can be qualified on exit,
!          FULL=.true. is returned. On entry, Diag is the local
!          diagonal. The global diagonal is synchronized if
!          Sync=.True. on entry.

use Cholesky, only: Cho_Real_Par, Diag_G, tMisc
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: Diag(*)
logical(kind=iwp), intent(in) :: Sync
integer(kind=iwp), intent(in) :: iShlAB, iSyMax, Mem
logical(kind=iwp), intent(out) :: Full
integer(kind=iwp) :: iLoc
real(kind=wp) :: c1, c2, w1, w2

call CWTime(c1,w1)

if (Cho_Real_Par) then

  ! Sync global diagonal (if requested).
  ! ------------------------------------

  if (Sync) then
    iLoc = 2
    call Cho_P_SyncDiag(Diag,iLoc)
  end if

  ! Swap local and global index arrays and use the original serial
  ! routines for qualifying diagonals.
  ! --------------------------------------------------------------
  !-TODO/FIXME: the memory constraint is then based on the global
  !             dimension, thus ensuring that all nodes find the same
  !             qualified (and this is absolutely essential!). However,
  !             one might find a more clever way of ensuring this....

  call Cho_P_IndxSwp()
  call Cho_Qualify(Diag_G,iShlAB,iSyMax,Mem,Full)
  call Cho_P_IndxSwp()

else

  call Cho_Qualify(Diag,iShlAB,iSyMax,Mem,Full)

end if

call CWTime(c2,w2)
tMisc(1,1) = tMisc(1,1)+c2-c1
tMisc(2,1) = tMisc(2,1)+w2-w1

end subroutine Cho_P_Qualify
