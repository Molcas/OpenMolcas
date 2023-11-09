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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************
! *******************************************
! ** Front-end routines for INPUT and MAIN **
! *******************************************

subroutine casvb(ireturn)

use Definitions, only: iwp, u5

implicit none
integer(kind=iwp), intent(out) :: ireturn

call cvbinp_cvb(0,u5)
call cvbmn_cvb(0)
call make_close_cvb(1)
ireturn = 0

return

end subroutine casvb
