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

subroutine Cho_P_SetRed(Diag,Sync)
!
! Purpose: set next reduced set after synchronizing the global
!          diagonal (if requested through flag Sync).
!          Global as well as local reduced sets are set.
!          Diag is the local diagonal, whereas Diag_G (defined in
!          Cholesky) points to the global diagonal.
!          Note that Diag is not referenced if Sync=.False.

use Cholesky, only: Cho_Real_Par, Diag_G
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: Diag(*)
logical(kind=iwp), intent(in) :: Sync
integer(kind=iwp) :: iLoc

if (Cho_Real_Par) then

  ! Sync global diagonal.
  ! ---------------------

  if (Sync) then
    iLoc = 2
    call Cho_P_SyncDiag(Diag,iLoc)
  end if

  ! Set next global reduced set. The original serial routines
  ! are used and so we must trick them by first swapping global
  ! and local index arrays (and then swap back, of course).
  ! -----------------------------------------------------------

  call Cho_P_IndxSwp()
  call Cho_SetRed(Diag_G)
  call Cho_P_IndxSwp()

  ! Set next local reduced set.
  ! ---------------------------

  call Cho_P_SetRed_L()

else

  call Cho_SetRed(Diag)

end if

end subroutine Cho_P_SetRed
