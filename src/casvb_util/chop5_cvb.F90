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

subroutine chop5_cvb()

use casvb_global, only: release

implicit real*8(a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

if (release(5)) call mfreer_cvb(ls(1))
release(5) = .true.
release(6) = .false.

ls(1) = mstackr_cvb(nsyme*norb*norb)
ls(2) = mstacki_cvb(ndimrel)
ls(3) = mstacki_cvb(norb)
ls(4) = mstackr_cvb(norb*norb*norb)
ls(5) = mstacki_cvb(2*(norb-1))
ls(6) = mstackr_cvb(min(norb-1,norbrel)*norb*norb)
ls(8) = mstacki_cvb(norb)
ls(9) = mstacki_cvb(nvb)
ls(10) = mstacki_cvb(nzrvb)
ls(11) = mstacki_cvb(2*nort)
ls(12) = mstacki_cvb(2*ndrot)
ls(13) = mstacki_cvb(nsyme)
if (.not. orbfr_is_unit) then
  ls(14) = mstackr_cvb(nprorb*nprorb)
else
  ls(14) = mstackr_cvb(0)
end if
if (iconstruc == 2) then
  ls(15) = mstackr_cvb(nvb*nvb)
else
  ls(15) = mstackr_cvb(0)
end if
ls(16) = mstacki_cvb(norb*nzeta)

return

end subroutine chop5_cvb
