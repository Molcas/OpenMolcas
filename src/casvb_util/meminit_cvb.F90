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

! -- Initialization of casvb memory manager ---
subroutine meminit_cvb()

use casvb_global, only: memdebug, nfieldm
use Definitions, only: u6

implicit none

memdebug = .false.
nfieldm = 0
call setmem('trace=off')
call setmem('clear=off')
if (memdebug) then
  write(u6,*) ' Casvb memory handler initialized.'
  write(u6,*) ' Memory offsets : integer= ',0,' real= ',0
  write(u6,*) ' No. of fields in use :',nfieldm
end if

return

end subroutine meminit_cvb
