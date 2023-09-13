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

subroutine mfreer_cvb(ipoint)
! Memory allocator. Releases pointer.

implicit real*8(a-h,o-z)
#include "memman_cvb.fh"

if (memdebug) write(6,*) '     Enter mfreer: pointer :',ipoint
! Check if allocated using mstack:
do ifield=1,nfield
  if (iaddr(ifield) == ipoint) then
    do jfield=ifield,nfield
      ipoint_g = iaddr(jfield)-ioff_r
      if (memdebug) write(6,*) '     Release pointer :',iaddr(jfield)
      call getmem('casvb','FREE','REAL',ipoint_g,nword)
    end do
    nfield = ifield-1
    return
  end if
end do
! Allocated through mheap:
ipoint_g = ipoint-ioff_r
call getmem('casvb','FREE','REAL',ipoint_g,nword)

return

end subroutine mfreer_cvb
