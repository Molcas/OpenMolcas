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

subroutine real_cvb(arr,nmax,nread,ifc)

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use casvb_global, only: inputmode
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nmax, ifc
real(kind=wp), intent(inout) :: arr(nmax)
integer(kind=iwp), intent(inout) :: nread
logical(kind=iwp) :: done

call real_cvb_internal(arr)

return

! This is to allow type punning without an explicit interface
contains

#include "intent.fh"

subroutine real_cvb_internal(arr)
  real(kind=wp), target, intent(_OUT_) :: arr(*)
  integer(kind=iwp), pointer :: iarr(:)
  integer(kind=iwp) :: i, ierr, ifcuse
  if (inputmode == 2) then
    call c_f_pointer(c_loc(arr(1)),iarr,[nread])
    call gethr_cvb(iarr,nread)
    nullify(iarr)
  else
    nread = 0
    if (nmax > 0) then
      ! Treat first field differently
      ifcuse = mod(ifc,4)
      if (ifcuse >= 2) ifcuse = 2
      call popfield_cvb(ifcuse)
      call rdreal_cvb(arr(1),ierr)
      done = .false.
      if (ierr <= 0) then
        nread = nread+1

        ifcuse = mod(ifc,2)
        done = .true.
        do i=2,nmax
          call popfield_cvb(ifcuse)
          call rdreal_cvb(arr(i),ierr)
          if (ierr > 0) then
            done = .false.
            exit
          end if
          nread = nread+1
        end do
      end if
      if (.not. done) then
        ! Crash if invalid field and IFC +4:
        if ((ierr == 4) .and. (ifc >= 4)) then
          write(u6,*) ' Invalid field found while reading real!'
          call abend_cvb()
        end if
        call pushfield_cvb()
      end if
    end if
    if (inputmode == 1) then
      call c_f_pointer(c_loc(arr(1)),iarr,[nread])
      call sethr_cvb(iarr,nread)
      nullify(iarr)
    end if
  end if
end subroutine real_cvb_internal

end subroutine real_cvb
