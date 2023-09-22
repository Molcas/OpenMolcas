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

subroutine mrealloci_cvb(ipoint,nword)

use casvb_global, only: memdebug
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: ipoint, nword
#include "idbl_cvb.fh"
integer(kind=iwp) :: iraddr, nwordr

iraddr = (ipoint-1)/idbl+1
nwordr = (nword+idbl-1)/idbl
call mreallocr_cvb(iraddr,nwordr)
ipoint = (iraddr-1)*idbl+1
if (memdebug) write(u6,*) '   mrealloci : nword & pointer :',nword,ipoint

return

end subroutine mrealloci_cvb
