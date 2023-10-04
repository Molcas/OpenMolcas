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

use casvb_global, only: civb1, civb2, civb3, civb4, civb5, civb6, civb7, civb8, civbvec, civbvecs, ndres_ok, release
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: iwp

implicit none
#include "main_cvb.fh"
#include "casinfo_cvb.fh"
integer(kind=iwp) :: iv

if (release(4)) then
  call mma_deallocate(civbvecs)
  nullify(civbvec)
  nullify(civb1)
  nullify(civb2)
  nullify(civb3)
  nullify(civb4)
  nullify(civb5)
  nullify(civb6)
  nullify(civb7)
  nullify(civb8)
end if
release(4) = .true.
release(5) = .false.

!FIXME: This deallocation should not be needed
if (allocated(civbvecs)) call mma_deallocate(civbvecs)
! CIVECP and CIVBH share memory
call mma_allocate(civbvecs,[0,ndres],[1,nv],label='civbvecs')
if (nv >= 1) then
  civb1 => civbvecs(1:,1)
  if (nv >= 2) then
    civb2 => civbvecs(1:,2)
  else
    civb2 => civbvecs(1:,1)
  end if
  if (nv >= 3) then
    civb3 => civbvecs(1:,3)
  else
    civb3 => civbvecs(1:,1)
  end if
  if (nv >= 4) then
    civb4 => civbvecs(1:,4)
  else
    civb4 => civbvecs(1:,1)
  end if
  if (nv >= 5) then
    civb5 => civbvecs(1:,5)
  else
    civb5 => civbvecs(1:,1)
  end if
  if (nv >= 6) then
    civb6 => civbvecs(1:,6)
  else
    civb6 => civbvecs(1:,1)
  end if
  if (nv >= 7) then
    civb7 => civbvecs(1:,7)
  else
    civb7 => civbvecs(1:,1)
  end if
  if (nv >= 8) then
    civb8 => civbvecs(1:,8)
  else
    civb8 => civbvecs(1:,1)
  end if
end if
civbvecs(0,:) = Zero
! Fix to put in "objects":
civbvec => civbvecs(2:,:)
do iv=1,nv
  call creatci_cvb(iv)
  if (.not. ndres_ok) call setcnt_cvb(civbvecs(1,iv),0)
end do
!-- ins
if (((ifinish == 1) .or. (ifinish == 2)) .and. (.not. lciweights)) then
  if (((.not. lcalcevb) .or. lcalccivbs) .and. ((.not. lcalcsvb) .or. lcalccivbs)) then
    civb3 => civbvecs(1:,1)
    civb4 => civbvecs(1:,1)
  else if ((.not. lcalcevb) .or. lcalccivbs) then
    civb4 => civbvecs(1:,3)
  else if ((.not. lcalcsvb) .or. lcalccivbs) then
    civb3 => civbvecs(1:,1)
    civb4 => civbvecs(1:,3)
  end if
end if
if (.not. memplenty) then
  select case (nv)
    case (2)
      civb3 => civbvecs(1:,1)
    case (3)
      civb4 => civbvecs(1:,1)
    case (4)
      civb5 => civbvecs(1:,1)
    case (5)
      civb6 => civbvecs(1:,1)
    case (6)
      civb7 => civbvecs(1:,1)
    case (7)
      civb8 => civbvecs(1:,1)
  end select
end if
if ((nv == 3) .or. (nv == 4) .or. (nv == 5)) then
  civb6 => civb2
  civb7 => civb3
  civb8 => civb4
end if

return

end subroutine chop4_cvb
