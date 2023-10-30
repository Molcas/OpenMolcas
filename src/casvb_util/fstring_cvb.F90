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

subroutine fstring_cvb(strings,nstring,istring,ncmp,ifc)

use casvb_global, only: inputmode
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nstring, ncmp, ifc
character(len=*), intent(in) :: strings(nstring)
integer(kind=iwp), intent(inout) :: istring
integer(kind=iwp) :: ierr
logical(kind=iwp) :: done
character(len=8) :: string
logical(kind=iwp), parameter :: debug = .false.

if (inputmode == 2) then
  call gethfs_cvb(istring)
  if (debug) then
    write(u6,*) ' fstring :',istring
    if (istring /= 0) write(u6,*) ' fstring :',strings(istring)
  end if
  return
end if
call popfield_cvb(ifc)
call rdstring_cvb(string,ierr)
done = .false.
do istring=1,nstring
  if (string(1:ncmp) == strings(istring)(1:ncmp)) then
    ! For keywords starting by END we check more letters. This
    ! implementation is a bit ungainly, but simpler to code:
    if (string(1:3) == 'END') then
      if (string(4:ncmp+3) /= strings(istring)(4:ncmp+3)) cycle
    end if
    done = .true.
    exit
  end if
end do
if (.not. done) then
  istring = 0
  call pushfield_cvb()
end if
if (inputmode == 1) call seth_cvb([istring],1)
if (debug) then
  write(u6,*) ' fstring :',istring
  if (istring /= 0) write(u6,*) ' fstring :',strings(istring)
end if

return

end subroutine fstring_cvb
