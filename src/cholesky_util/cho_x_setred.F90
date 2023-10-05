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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************
!  Cho_X_SetRed
!
!> @brief
!>   Read and set index arrays for reduced set \p iRed at location \p iLoc
!> @author Thomas Bondo Pedersen
!>
!> @details
!> Reads information for reduced set \p iRed (= ``1``, ``2``, ..., \c MaxRed)
!> and sets up the index arrays
!>
!> - \p nnBstRT(iLoc)      &rarr; accessible via Cholesky
!> - \p nnBstR(:,iLoc)     &rarr; accessible via Cholesky
!> - \p iiBstR(:,iLoc)     &rarr; accessible via Cholesky
!> - \p nnBstRSh(:,:,iLoc) &rarr; accessible via Cholesky
!> - \p iiBstRSh(:,:,iLoc) &rarr; accessible via Cholesky
!> - \p IndRed(:,iLoc)     &rarr; accessible via Cholesky
!>
!> On succesful completion, \p irc = ``0`` is returned.
!> Note that the only allowed \p iLoc values are ``2`` and ``3``; any other
!> value gives rise to error code \p irc = ``1`` and nothing is done!
!> If \p iRed is out of bounds, \p irc = ``2`` is returned and nothing is done!
!>
!> @note
!> The Cholesky procedure must have been successfully initialized (by ::Cho_X_Init).
!>
!> @param[out] irc  return code
!> @param[in]  iLoc location in index arrays
!> @param[in]  iRed reduced set on disk
!***********************************************************************

subroutine Cho_X_SetRed(irc,iLoc,iRed)

use Cholesky, only: IndRed, MaxRed
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: iLoc, iRed
integer(kind=iwp) :: iab

if ((iLoc == 2) .or. (iLoc == 3)) then
  if ((iRed < 1) .or. (iRed > MaxRed)) then
    irc = 2
  else
    call Cho_GetRed(iRed,iLoc,.false.)
    call Cho_SetRedInd(iLoc)
    irc = 0
    if (iRed == 1) then ! set correct IndRed array
      do iab=1,size(IndRed,1)
        IndRed(iab,iLoc) = iab
      end do
    end if
  end if
else
  irc = 1
end if

end subroutine Cho_X_SetRed
