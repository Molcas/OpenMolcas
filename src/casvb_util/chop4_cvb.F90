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

subroutine chop4_cvb()

use casvb_global, only: ndres_ok, release
use Constants, only: Zero
use Definitions, only: iwp

implicit none
#include "main_cvb.fh"
#include "casinfo_cvb.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: iv
integer(kind=iwp), external :: mstackr_cvb

if (release(4)) call mfreer_cvb(lc(1)-2)
release(4) = .true.
release(5) = .false.

! CIVECP and CIVBH share memory --> LC(2)
do iv=1,nv
  lc(iv) = mstackr_cvb(ndres)
end do
do iv=1,nv
  work(lc(iv)) = Zero
  work(lc(iv)+1) = Zero
end do
do iv=1,nv
  lc(iv) = lc(iv)+2
end do
! Fix to put in "objects":
do iv=1,nv
  call creatci_cvb(iv,work(lc(iv)),lc(iv)+1,nint(work(lc(iv)-2)),work(lc(iv)-1))
  if (.not. ndres_ok) call setcnt_cvb(work(lc(iv)),0)
end do
!-- ins
if (((ifinish == 1) .or. (ifinish == 2)) .and. (.not. lciweights)) then
  if (((.not. lcalcevb) .or. lcalccivbs) .and. ((.not. lcalcsvb) .or. lcalccivbs)) then
    lc(3) = lc(1)
    lc(4) = lc(1)
  else if ((.not. lcalcevb) .or. lcalccivbs) then
    lc(4) = lc(3)
  else if ((.not. lcalcsvb) .or. lcalccivbs) then
    lc(4) = lc(3)
    lc(3) = lc(1)
  end if
end if
if (.not. memplenty) lc(nv+1) = lc(1)
if ((nv == 3) .or. (nv == 4) .or. (nv == 5)) then
  lc(6) = lc(2)
  lc(7) = lc(3)
  lc(8) = lc(4)
end if

return

end subroutine chop4_cvb
