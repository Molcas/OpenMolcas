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
! -- All memory allocated via real*8 routines ---

integer function mheapr_cvb(nword1)
! Memory allocator (heap). Returns pointer to NWORD real*8 words.

use casvb_global, only: ioff_r, memdebug

implicit real*8(a-h,o-z)

nword = nword1
if (memdebug) write(6,*) '     Enter mheapr: nword :',nword
if (nword < 0) then
  write(6,*) ' Error: attempting to allocate negative amount of memory.'
  write(6,*) ' nword=',nword
  call abend_cvb()
end if
call getmem('casvb','ALLO','REAL',ipoint_g,nword)
mheapr_cvb = ipoint_g+ioff_r
if (memdebug) write(6,*) '     mheapr: nword & pointer :',nword,mheapr_cvb

return

end function mheapr_cvb
