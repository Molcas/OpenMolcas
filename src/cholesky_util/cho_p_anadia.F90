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

subroutine Cho_P_AnaDia(Diag,Sync,Bin1,Step,NumBin,Full)
!
! Purpose: analyze global diagonal (histogram). Diag is the local
!          diagonal. If Sync=.True. the global diagonal is
!          synchronized before analysis.

use Cholesky, only: Cho_Real_Par, Diag_G
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: Diag(*), Bin1, Step
logical(kind=iwp), intent(in) :: Sync, Full
integer(kind=iwp) :: NumBin
integer(kind=iwp) :: iLoc

if (Cho_Real_Par) then

  ! If requested, sync diagonal.
  ! ----------------------------

  if (Sync) then
    iLoc = 2
    call Cho_P_SyncDiag(Diag,iLoc)
  end if

  ! Swap local and global index arrays and use original serial
  ! serial to perform analysis of the global diagonal.
  ! ----------------------------------------------------------

  call Cho_P_IndxSwp()
  call Cho_AnaDia(Diag_G,Bin1,Step,NumBin,Full)
  call Cho_P_IndxSwp()

else

  call Cho_AnaDia(Diag,Bin1,Step,NumBin,Full)

end if

end subroutine Cho_P_AnaDia
