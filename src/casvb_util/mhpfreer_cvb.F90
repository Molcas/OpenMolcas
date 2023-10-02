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

subroutine mhpfreer_cvb(ipoint)
! Memory allocator. Releases pointer.

use casvb_global, only: iaddrm, memdebug, nfieldm
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: ipoint
integer(kind=iwp) :: ifield, ipoint_g, jfield, nword

if (memdebug) write(u6,*) '     Enter mfreer: pointer :',ipoint
! Check if allocated using mstack:
do ifield=1,nfieldm
  if (iaddrm(ifield) == ipoint) then
    do jfield=ifield,nfieldm
      ipoint_g = iaddrm(jfield)
      if (memdebug) write(u6,*) '     Release pointer :',iaddrm(jfield)
      call getmem('casvb','FREE','REAL',ipoint_g,nword)
    end do
    nfieldm = ifield-1
    return
  end if
end do
! Allocated through mheap:
ipoint_g = ipoint
call getmem('casvb','FREE','REAL',ipoint_g,nword)

return

end subroutine mhpfreer_cvb
