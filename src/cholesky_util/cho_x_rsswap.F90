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
!  Cho_X_RSSwap
!
!> @brief
!>   Swap reduced set index arrays
!> @author Thomas Bondo Pedersen
!>
!> @details
!> Swap index arrays at locations \p iRS (``1``,``2``,``3``) and \p jRS (``1``,``2``,``3``).
!> If \p iRS and/or \p jRS are out of bounds, \p irc = ``1`` on exit and nothing
!> has been done!
!>
!> @warning
!> No special action is taken to redefine the \c IndRed array for first reduced set.
!>
!> @note
!> The CHolesky module must have been initialized.
!>
!> @param[out] irc return code
!> @param[in]  iRS location of reduced set
!> @param[in]  jRS location of reduced set
!***********************************************************************

subroutine Cho_X_RSSwap(irc,iRS,jRS)

use Cholesky, only: iiBstR, iiBstRSh, IndRed, nnBstR, nnBstRSh, nnBstRT, nnShl, nSym
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: iRS, jRS
integer(kind=iwp) :: N, iTemp

if ((iRS < 1) .or. (iRS > 3) .or. (jRS < 1) .or. (jRS > 3)) then
  irc = 1
else
  if (iRS /= jRS) then
    N = nSym*nnShl
    call iSwap(N,iiBstRSh(:,:,iRS),1,iiBstRSh(:,:,jRS),1)
    call iSwap(N,nnBstRSh(:,:,iRS),1,nnBstRSh(:,:,jRS),1)
    call iSwap(nSym,iiBstR(1,iRS),1,iiBstR(1,jRS),1)
    call iSwap(nSym,nnBstR(1,iRS),1,nnBstR(1,jRS),1)
    call iSwap(nnBstRT(1),IndRed(:,iRs),1,IndRed(:,jRs),1)
    iTemp = nnBstRT(iRS)
    nnBstRT(iRS) = nnBstRT(jRS)
    nnBstRT(jRS) = iTemp
  end if
  irc = 0
end if

end subroutine Cho_X_RSSwap
