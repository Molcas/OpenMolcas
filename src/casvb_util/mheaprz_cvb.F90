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

integer function mheaprz_cvb(nword)

implicit real*8(a-h,o-z)
#include "memman_cvb.fh"
#include "WrkSpc.fh"

if (memdebug) write(6,*) ' mheaprz :'
mheaprz_cvb = mheapr_cvb(nword)
call fzero(work(mheaprz_cvb),nword)

return

end function mheaprz_cvb
