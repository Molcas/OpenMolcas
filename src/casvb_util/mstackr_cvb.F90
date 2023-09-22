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

function mstackr_cvb(nword)
! Memory allocator (stack). Returns pointer to NWORD real*8 words.

use casvb_global, only: iaddrm, memdebug, mxfield, nfieldm
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: mstackr_cvb
integer(kind=iwp) :: nword
integer(kind=iwp), external :: mheapr_cvb

if (memdebug) write(u6,*) '     Enter mstackr: nword :',nword
mstackr_cvb = mheapr_cvb(nword)
nfieldm = nfieldm+1
if (nfieldm > mxfield) then
  write(u6,*) ' Too many field in mstackr :',nfieldm,mxfield
  call abend_cvb()
end if
iaddrm(nfieldm) = mstackr_cvb
if (memdebug) write(u6,*) '     mstackr: nword & pointer :',nword,mstackr_cvb

return

end function mstackr_cvb
