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

subroutine mreallocr_cvb(ipoint,nword)
! Memory allocator (heap). Reallocate pointer.

implicit real*8(a-h,o-z)
#include "memman_cvb.fh"
#include "WrkSpc.fh"
#include "files_cvb.fh"

if (memdebug) write(6,*) '     Enter mreallocr: nword & pointer :',nword,ipoint

ipoint_g = ipoint-ioff_r

!call getmem('casvb','CHAN','REAL',ipoint_g,nword)

! Read and write data -- not efficient but safe and simple:
call getmem('casvb','LENG','REAL',ipoint_g,nword_old)
nword_move = min(nword,nword_old)
call wrr_cvb(work(ipoint),nword_move,recn_tmp06,0)
call mfreer_cvb(ipoint)
ipoint = mheapr_cvb(nword)
call rdr_cvb(work(ipoint),nword_move,recn_tmp06,0)

if (memdebug) write(6,*) '     mreallocr : nword & pointer :',nword,ipoint

return

end subroutine mreallocr_cvb
