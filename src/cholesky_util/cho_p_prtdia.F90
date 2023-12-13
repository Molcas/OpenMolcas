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

subroutine Cho_P_PrtDia(Diag,Sync,iSyLst,nSyLst,iLoc)
!
! Purpose: print global diagonal. Diag is the local diagonal and
!          if Sync=.True. the global diagonal is
!          synchronized before printing. Array iSyLst(nSyLst)
!          specifies which symmetry blocks to print, and iLoc points
!          to the memory location of the reduced set index arrays to
!          use for printing (and synchronizing, if requested).

use Cholesky, only: Cho_Real_Par, Diag_G
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: Diag(*)
logical(kind=iwp), intent(in) :: Sync
integer(kind=iwp), intent(in) :: nSyLst, iSyLst(nSyLst), iLoc

if (Cho_Real_Par) then

  ! Sync diagonal if requested.
  ! ---------------------------

  if (Sync) call Cho_P_SyncDiag(Diag,iLoc)

  ! Swap local and global index arrays and use original serial routine
  ! to print diagonal.
  ! ------------------------------------------------------------------

  call Cho_P_IndxSwp()
  call Cho_PrtDia(Diag_G,iSyLst,nSyLst,iLoc)
  call Cho_P_IndxSwp()

else

  call Cho_PrtDia(Diag,iSyLst,nSyLst,iLoc)

end if

end subroutine Cho_P_PrtDia
