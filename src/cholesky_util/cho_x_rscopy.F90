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
!  Cho_X_RSCopy
!
!> @brief
!>   Copy reduced set index arrays
!> @author Thomas Bondo Pedersen
!>
!> @details
!> Copy index arrays from location \p iRS (``1``,``2``,``3``) to location
!> \p jRS (``1``,``2``,``3``). Special handling of the case \p iRS = ``1``:
!> ``IndRed(i,jRS) = i`` rather than ``IndRed(i,jRS)=IndRed(i,iRS)``.
!> If \p iRS and/or \p jRS are out of bounds, \p irc = ``1`` on exit and nothing
!> has been done!
!>
!> @note
!> The Cholesky module must have been initialized.
!>
!> @param[out] irc return code
!> @param[in]  iRS location of reference reduced set
!> @param[in]  jRS location of target reduced set
!***********************************************************************

subroutine Cho_X_RSCopy(irc,iRS,jRS)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: iRS, jRS

if ((iRS < 1) .or. (iRS > 3) .or. (jRS < 1) .or. (jRS > 3)) then
  irc = 1
else
  call Cho_RSCopy(iRS,jRS)
  irc = 0
end if

end subroutine Cho_X_RSCopy
