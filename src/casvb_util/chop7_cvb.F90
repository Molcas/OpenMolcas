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

subroutine chop7_cvb()

implicit real*8(a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
#include "rls_cvb.fh"
#include "change7.fh"

if (release(7)) call mfreer_cvb(lw(1))
release(7) = .true.
release(8) = .false.

if (icase == 1) then
  lw(1) = mstackr_cvb(norb*norb)
  lw(2) = mstackr_cvb(norb*norb)
  lw(3) = mstackr_cvb(norb*norb)

  lw(4) = mstackr_cvb(norb*norb+ihlf_cvb(norb+2*norb*norb))
  lw(5) = mstackr_cvb(norb*norb+ihlf_cvb(norb+2*norb*norb))
  lw(6) = mstackr_cvb(norb*norb+ihlf_cvb(norb+2*norb*norb))

  lw(7) = mstackr_cvb(nvb)
  lw(8) = mstackr_cvb(nvb)
  lw(9) = mstackr_cvb(ndetvb)
  lw(10) = mstackr_cvb(ndetvb)
  lw(11) = mstackr_cvb(ndetvb)
else if (icase == 2) then
  lw(1) = mstackr_cvb(norb*norb)
  lw(2) = mstackr_cvb(norb*norb)
  lw(3) = mstackr_cvb(norb*norb)

  lw(4) = mstackr_cvb(norb*norb+ihlf_cvb(norb+2*norb*norb))
  lw(5) = mstackr_cvb(norb*norb+ihlf_cvb(norb+2*norb*norb))
  lw(6) = mstackr_cvb(norb*norb+ihlf_cvb(norb+2*norb*norb))

  lw(7) = mstackr_cvb(nvb)
  lw(8) = mstackr_cvb(nvb)
  lw(9) = mstackr_cvb(ndetvb)
  lw(10) = mstackr_cvb(ndetvb)
  lw(11) = mstackr_cvb(ndetvb)
else if (icase == 3) then
  lw(1) = mstackr_cvb(norb*norb)
  lw(2) = mstackr_cvb(norb*norb)
  lw(3) = mstackr_cvb(norb*norb)

  lw(4) = mstackr_cvb(norb*norb+ihlf_cvb(norb+2*norb*norb))
  lw(5) = mstackr_cvb(norb*norb+ihlf_cvb(norb+2*norb*norb))
  lw(6) = mstackr_cvb(norb*norb+ihlf_cvb(norb+2*norb*norb))

  lw(7) = mstackr_cvb(nvb)
  lw(8) = mstackr_cvb(nvb)
  lw(9) = mstackr_cvb(ndetvb)
  lw(10) = mstackr_cvb(ndetvb)
  lw(11) = mstackr_cvb(ndetvb)
end if
if (anyslater) then
  nvb_alloc = ndetvb
else
  nvb_alloc = ndetvb
end if
lw(12) = mstackr_cvb(norb*norb)
lw(13) = mstackr_cvb(nvb_alloc)

lv(5) = mstackr_cvb(ndetvb)

ibasemx = max(ibasemx,mstackr_cvb(0))

return

end subroutine chop7_cvb
