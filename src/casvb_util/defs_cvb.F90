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

subroutine defs_cvb()

use casvb_global, only: iunset

implicit real*8(a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
logical :: ifploc = .false.

! Default settings:
strtvb = zero
savvb = zero
savvbci = zero
kbasis = 1
mxiter = iunset
icrit = iunset
imethod = iunset
isaddle = iunset
initial = -1
opposite = .false.
projcas = .false.
projsym = .false.
npcf = iunset
ishstruc = iunset
! +1=CHIRGWIN +2=LOWDIN +4=INVERSE
ivbweights = iunset
iciweights = iunset
opposite = .false.
projcas = .false.
projsym = .false.
sij = .false.
anyslater = .false.
service = .false.
do i=1,10
  ip(i) = 1
end do

call tunedefs_cvb()
if (variat) then
  ploc = ifploc
else
  ploc = .false.
end if

! Counters
norbrel = 0
ndimrel = 0
nfxorb = 0
nfxvb = 0
nzrvb = 0
nort = 0
lfxvb = 0
lzrvb = 0

return

end subroutine defs_cvb
