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
! -- Zeroing routines - just front-ends ---

function mstackrz_cvb(nword)

use casvb_global, only: memdebug
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: mstackrz_cvb
integer(kind=iwp) :: nword
#include "WrkSpc.fh"
integer(kind=iwp), external :: mstackr_cvb

if (memdebug) write(u6,*) ' mstackrz :'
mstackrz_cvb = mstackr_cvb(nword)
call fzero(work(mstackrz_cvb),nword)

return

end function mstackrz_cvb
