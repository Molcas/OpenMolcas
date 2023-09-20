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

subroutine o7a_cvb(nparm)

use casvb_global, only: have_solved_it

implicit real*8(a-h,o-z)
save one
data one/1d0/

call ddnewopt_cvb()
have_solved_it = .false.
call ddguess_cvb([one],1,0)

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(nparm)

end subroutine o7a_cvb
