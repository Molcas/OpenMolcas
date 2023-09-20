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

use casvb_global, only: ioff_i, ioff_r, memdebug, nfieldm

implicit real*8(a-h,o-z)

memdebug = .false.
nfieldm = 0
ioff_r = 0
ioff_i = 0
call setmem('trace=off')
call setmem('clear=off')
if (memdebug) then
  write(6,*) ' Casvb memory handler initialized.'
  write(6,*) ' Memory offsets : integer= ',ioff_i,' real= ',ioff_r
  write(6,*) ' No. of fields in use :',nfieldm
end if

return

end subroutine meminit_cvb
