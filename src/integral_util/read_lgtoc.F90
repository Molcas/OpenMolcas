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
! Copyright (C) 1992, Roland Lindh                                     *
!***********************************************************************

subroutine read_lgtoc(lgtoc,gtoc,n)
!*********** columbus interface ****************************************
!read table of contents for gamma file

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n, lgtoc
real(kind=wp), intent(out) :: gtoc(n)

read(lgtoc) gtoc

return

end subroutine read_lgtoc
