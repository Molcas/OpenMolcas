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

subroutine casinfodef_cvb()

use casvb_global, only: inputmode, iorclos_d, iorcore_d, iorocc_d, mxorb_cvb, noe, nstsym_d, recn_jobiph, recn_jobold, &
                        recn_oneint, recn_vbwfn, strtci, strtint, strtmo, strtvb, savvb, savvbci, variat

implicit none

! Counters
nstsym_d = 0

noe = 2*mxorb_cvb
if (.not. variat) then
  strtci = recn_jobold
else
  strtci = recn_jobiph
end if
strtmo = recn_jobiph
strtint = recn_oneint
strtvb = recn_vbwfn
savvb = recn_vbwfn
savvbci = recn_jobiph

if (inputmode == 2) then
  iorcore_d(:) = -1
  iorclos_d(:) = -1
  iorocc_d(:) = -1
end if

return

end subroutine casinfodef_cvb
