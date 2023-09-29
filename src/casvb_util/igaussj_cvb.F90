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

subroutine igaussj_cvb(orbs,igjorb)

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
real(kind=wp) :: orbs(norb,norb)
integer(kind=iwp) :: igjorb(*)
integer(kind=iwp) :: i, ioff
integer(kind=iwp), allocatable :: lrow(:)
real(kind=wp), allocatable :: a(:,:)
integer(kind=iwp), external :: idbl_cvb

call mma_allocate(a,norb,norb,label='a')
call mma_allocate(lrow,norb,label='lrow')
call fmove_cvb(orbs,a,norb*norb)
ioff = idbl_cvb(norb*norb)
call igaussj_cvb_internal(igjorb)
call imove_cvb(igjorb(1+ioff),lrow,norb)
do i=1,norb
  igjorb(lrow(i)+ioff) = i
end do
call mma_deallocate(a)
call mma_deallocate(lrow)

return

! This is to allow type punning without an explicit interface
contains

subroutine igaussj_cvb_internal(igjorb)
  integer(kind=iwp), target :: igjorb(*)
  real(kind=wp), pointer :: gjorb(:)
  call c_f_pointer(c_loc(igjorb(1)),gjorb,[norb*norb])
  call gaussj2_cvb(a,lrow,igjorb(1+ioff),igjorb(1+norb+ioff),gjorb,norb)
  nullify(gjorb)
end subroutine igaussj_cvb_internal

end subroutine igaussj_cvb
