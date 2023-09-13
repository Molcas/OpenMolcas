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

implicit real*8(a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
#include "WrkSpc.fh"
dimension orbs(norb,norb), igjorb(*)

call igaussj_cvb_internal(igjorb)

! This is to allow type punning without an explicit interface
contains

subroutine igaussj_cvb_internal(igjorb)
  use iso_c_binding
  integer, target :: igjorb(*)
  real*8, pointer :: gjorb(:)
  k1 = mstackr_cvb(norb*norb)
  k2 = mstacki_cvb(norb)
  k3 = mstacki_cvb(norb)
  k4 = mstacki_cvb(norb)
  call fmove_cvb(orbs,work(k1),norb*norb)
  ioff = idbl_cvb(norb*norb)
  call c_f_pointer(c_loc(igjorb(1)),gjorb,[norb*norb])
  call gaussj2_cvb(work(k1),iwork(k2),iwork(k3),iwork(k4),igjorb(1+ioff),igjorb(1+norb+ioff),gjorb,norb)
  nullify(gjorb)
  call imove_cvb(igjorb(1+ioff),iwork(k2),norb)
  do i=1,norb
    igjorb(iwork(i+k2-1)+ioff) = i
  end do
  call mfreer_cvb(k1)
  return
end subroutine igaussj_cvb_internal

end subroutine igaussj_cvb
