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

use casvb_global, only: inputmode

implicit real*8(a-h,o-z)
dimension arr(nmax)
logical done

call real_cvb_internal(arr)

! This is to allow type punning without an explicit interface
contains

subroutine real_cvb_internal(arr)
  use iso_c_binding
  real*8, target :: arr(*)
  integer, pointer :: iarr(:)
  if (inputmode == 2) then
    call c_f_pointer(c_loc(arr(1)),iarr,[nread])
    call gethr_cvb(iarr,nread)
    nullify(iarr)
    return
  end if
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
        write(6,*) ' Invalid field found while reading real!'
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
  return
end subroutine real_cvb_internal

end subroutine real_cvb
