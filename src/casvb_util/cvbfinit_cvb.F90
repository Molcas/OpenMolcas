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

subroutine cvbfinit_cvb()

use casvb_global, only: is_set
use Definitions, only: iwp, RtoI

implicit none
#include "Molcas.fh"
#include "main_cvb.fh"
#include "idbl_cvb.fh"
#include "print_cvb.fh"
integer(kind=iwp), parameter :: iset = 1

mxaobf = maxbfn
iprec = 8
iwidth = 110
call formats_cvb()
idbl = RtoI
call setmem('trace=off')
call setmem('clear=off')
if (is_set /= iset) then
  ! Initializations below are only carried out once:
  call io_init_cvb()
  call main_bdata_cvb()
end if
is_set = iset

return

end subroutine cvbfinit_cvb
