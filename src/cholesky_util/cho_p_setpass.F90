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

subroutine Cho_P_SetPass(Diag,Sync,DiaSh,iSySh,iLoc,Conv,nPotSh)
!
! Purpose: check convergence and, if not converged, set up next
!          next integral pass
!          Diag is the local diagonal; the global diagonal is
!          synchronized if Sync=.True. Note that DiaSh and iSySh
!          must be allocated with dimension nnShl_G (global number
!          of shell pairs in 1st reduced set). iLoc is the location
!          to use in the reduced set index arrays. On exit,
!          Conv=.True. if converged and nPotSh is the number of
!          shell pairs whose max. diagonal element is larger than
!          the decomposition threshold.

use Cholesky, only: Cho_Real_Par, Diag_G
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: Diag(*)
logical(kind=iwp), intent(in) :: Sync
real(kind=wp), intent(_OUT_) :: DiaSh(*)
integer(kind=iwp), intent(_OUT_) :: iSySh(*)
integer(kind=iwp), intent(in) :: iLoc
logical(kind=iwp), intent(out) :: Conv
integer(kind=iwp), intent(out) :: nPotSh

if (Cho_Real_Par) then

  ! Sync diagonal if requested.
  ! ---------------------------

  if (Sync) call Cho_P_SyncDiag(Diag,iLoc)

  ! Swap local and global index arrays and set next integral pass
  ! original serial routine.
  ! -------------------------------------------------------------

  call Cho_P_IndxSwp()
  call Cho_SetPass(Diag_G,DiaSh,iSySh,iLoc,Conv,nPotSh)
  call Cho_P_IndxSwp()

else

  call Cho_SetPass(Diag,DiaSh,iSySh,iLoc,Conv,nPotSh)

end if

end subroutine Cho_P_SetPass
