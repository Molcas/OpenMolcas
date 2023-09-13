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

implicit real*8(a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
#include "casinfo_cvb.fh"
#include "inpmod_cvb.fh"

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
  do i=1,mxstsy_ci
    iorcore_d(i) = -1
    iorclos_d(i) = -1
    iorocc_d(i) = -1
  end do
end if

return

end subroutine casinfodef_cvb
