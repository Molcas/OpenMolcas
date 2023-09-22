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

function mheapr_cvb(nword1)
! Memory allocator (heap). Returns pointer to NWORD real*8 words.

use casvb_global, only: ioff_r, memdebug
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: mheapr_cvb
integer(kind=iwp) :: nword1
integer(kind=iwp) :: ipoint_g, nword

nword = nword1
if (memdebug) write(u6,*) '     Enter mheapr: nword :',nword
if (nword < 0) then
  write(u6,*) ' Error: attempting to allocate negative amount of memory.'
  write(u6,*) ' nword=',nword
  call abend_cvb()
end if
call getmem('casvb','ALLO','REAL',ipoint_g,nword)
mheapr_cvb = ipoint_g+ioff_r
if (memdebug) write(u6,*) '     mheapr: nword & pointer :',nword,mheapr_cvb

return

end function mheapr_cvb
