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
! -- Integer routines - just front-ends for real*8 ---

integer function mstacki_cvb(nword)

implicit real*8(a-h,o-z)
#include "idbl_cvb.fh"
#include "memman_cvb.fh"

if (memdebug) write(6,*) '   Enter mstacki: nword :',nword
nwordr = (nword+idbl-1)/idbl
iraddr = mstackr_cvb(nwordr)
mstacki_cvb = (iraddr-1)*idbl+1
if (memdebug) write(6,*) '   mstacki: nword & pointer :',nword,mstacki_cvb

return

end function mstacki_cvb
