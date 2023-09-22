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

!***********************************************************************
!*                                                                     *
!*  GAUSSJ    := Define sequence of simple updates from orb transf.    *
!*                                                                     *
!***********************************************************************
subroutine gaussj_cvb(orbs,gjorb)

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
real(kind=wp) :: orbs(norb,norb), gjorb(*)

call gaussj_cvb_internal(gjorb)

! This is to allow type punning without an explicit interface
contains

subroutine gaussj_cvb_internal(gjorb)
  real(kind=wp), target :: gjorb(*)
  integer(kind=iwp), pointer :: igjorb(:)
  call c_f_pointer(c_loc(gjorb(1)),igjorb,[1])
  call igaussj_cvb(orbs,igjorb)
  nullify(igjorb)
end subroutine gaussj_cvb_internal

end subroutine gaussj_cvb
