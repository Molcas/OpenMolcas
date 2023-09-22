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
!************************************
!** Routines involving CI and ORBS **
!************************************

subroutine iapplyt_cvb(cvec,igjorb)

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: cvec(*)
integer(kind=iwp) :: igjorb(*)
#include "main_cvb.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: ioff, ivec
integer(kind=iwp), external :: idbl_cvb

ivec = nint(cvec(1))
n_applyt = n_applyt+1
ioff = idbl_cvb(norb*norb)
if (iform_ci(ivec) == 0) then
  call permci_cvb(cvec,igjorb(1+ioff))
  call iapplyt_cvb_internal(igjorb)
else
  write(u6,*) ' Unsupported format in APPLYT :',iform_ci(ivec)
  call abend_cvb()
end if
call setcnt2_cvb(ivec,0)

return

! This is to allow type punning without an explicit interface
contains

subroutine iapplyt_cvb_internal(igjorb)
  integer(kind=iwp), target :: igjorb(*)
  real(kind=wp), pointer :: gjorb(:)
  call c_f_pointer(c_loc(igjorb(1)),gjorb,[1])
  call applyt2_cvb(work(iaddr_ci(ivec)),gjorb,igjorb(1+norb+ioff),iwork(ll(1)),iwork(ll(2)),iwork(ll(5)),iwork(ll(6)),work(ll(9)), &
                   work(ll(10)))
  nullify(gjorb)
end subroutine iapplyt_cvb_internal

end subroutine iapplyt_cvb
