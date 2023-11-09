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

subroutine string_cvb(arr,nmax,nread,ifc)

use casvb_global, only: inputmode
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nmax, ifc
character(len=*), intent(inout) :: arr(nmax)
integer(kind=iwp), intent(inout) :: nread
integer(kind=iwp) :: i, ierr, ifcuse
character(len=100) :: string
logical(kind=iwp) :: done

if (inputmode == 2) then
  call geths_cvb(arr,nread)
  return
end if
nread = 0

if (nmax > 0) then
  ! Treat first field differently
  ifcuse = mod(ifc,4)
  if (ifcuse >= 2) ifcuse = 2
  call popfield_cvb(ifcuse)
  call rdstring_cvb(string,ierr)
  done = .false.
  if (ierr <= 0) then
    arr(1) = string
    nread = nread+1

    ifcuse = mod(ifc,2)
    done = .true.
    do i=2,nmax
      call popfield_cvb(ifcuse)
      call rdstring_cvb(string,ierr)
      if (ierr > 0) then
        done = .false.
        exit
      end if
      arr(i) = string
      nread = nread+1
    end do
  end if
  if (.not. done) call pushfield_cvb()
end if
if (inputmode == 1) call seths_cvb(arr,nread)

return

end subroutine string_cvb
