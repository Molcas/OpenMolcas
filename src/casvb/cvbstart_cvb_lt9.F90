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

subroutine cvbstart_cvb_lt9(icode)

use casvb_global, only: endvar, iprm, iunset, lstprm, recinp, recinp_old, variat
use Constants, only: Zero
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: icode

variat = (mod(icode,10) /= 0)
endvar = (mod(icode,10) == 2)
recinp = Zero
recinp_old = Zero
lstprm(:) = iunset
iprm = 0
if (.not. variat) call casinfo1_cvb()

return

end subroutine cvbstart_cvb_lt9
