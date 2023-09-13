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

integer function mstackr_cvb(nword)
! Memory allocator (stack). Returns pointer to NWORD real*8 words.
implicit real*8(a-h,o-z)
#include "memman_cvb.fh"

if (memdebug) write(6,*) '     Enter mstackr: nword :',nword
mstackr_cvb = mheapr_cvb(nword)
nfield = nfield+1
if (nfield > mxfield) then
  write(6,*) ' Too many field in mstackr :',nfield,mxfield
  call abend_cvb()
end if
iaddr(nfield) = mstackr_cvb
if (memdebug) write(6,*) '     mstackr: nword & pointer :',nword,mstackr_cvb

return

end function mstackr_cvb
