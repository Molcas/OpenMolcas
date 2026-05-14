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

subroutine Cho_P_MaxDX(Diag,Sync,Dmax)
!
! Purpose: get max. diagonal elements in each sym. block,
!          qualified diagonals excluded.

use Cholesky, only: Cho_Real_Par, Diag_G
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(inout) :: Diag(*)
real(kind=wp), intent(_OUT_) :: Dmax(*)
logical(kind=iwp), intent(in) :: Sync
integer(kind=iwp) :: iLoc

if (Cho_Real_Par) then
  if (Sync) then
    iLoc = 2
    call Cho_P_SyncDiag(Diag,iLoc)
  end if
  call Cho_P_IndxSwp()
  call Cho_MaxDX(Diag_G,Dmax)
  call Cho_P_IndxSwp()
else
  call Cho_MaxDX(Diag,Dmax)
end if

end subroutine Cho_P_MaxDX
