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

subroutine cvbmn_cvb(icode)

use casvb_global, only: esym, n_iter
use Definitions, only: iwp
use lucia_interface, only: lucia_util

implicit none
integer(kind=iwp), intent(in) :: icode

! ICODE=0 standard casvb calculation
! ICODE=1 variational calculation
! ICODE=2 end of variational calculation (print summary)

call cvbstart_cvb_lt9(icode)
call main_cvb()
call setretvals_cvb(esym,n_iter)
call Lucia_Util('CLOSE')

end subroutine cvbmn_cvb
