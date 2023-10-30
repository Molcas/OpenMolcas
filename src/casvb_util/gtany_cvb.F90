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

subroutine gtany_cvb(string,intx,realx,ic,ifield,ierr)

use casvb_global, only: iline, ilv, lenline, line
use Definitions, only: wp, iwp, u6

implicit none
character(len=*), intent(inout) :: string
integer(kind=iwp), intent(inout) :: intx, ierr
real(kind=wp), intent(inout) :: realx
integer(kind=iwp), intent(in) :: ic, ifield
integer(kind=iwp) :: ich, iempty, ifirst, istatus, jch, jfield, jline
logical(kind=iwp) :: done
integer(kind=iwp), parameter :: nempty = 1
logical(kind=iwp), parameter :: debug = .false.
character(len=*), parameter :: empty(nempty) = ['--']
logical(kind=iwp), external :: isitanint_cvb, isitareal_cvb

if (ic > 1) ierr = 0
jfield = 1
jline = 1
do ich=1,lenline
  if (ilv(ich) == 1) jline = jline+1
  if ((jline == iline) .and. (ilv(ich) == 2)) jfield = jfield+1
  if ((jline == iline) .and. (jfield == ifield)) then
    jch = ich
    if (ich == 1) jch = 0
    do
      jch = jch+1
      if ((ilv(jch) /= 0) .or. (jch == lenline+1)) exit
    end do
    ifirst = ich+1
    if (ich == 1) ifirst = 1
    ! Special character strings to signify empty field?
    done = .false.
    do iempty=1,nempty
      if (line(ifirst:jch-1) == trim(empty(iempty))) then
        done = .true.
        exit
      end if
    end do
    if (done) then
      ! "Empty" field:
      if (ic == 1) then
        string = ' '
      else
        ierr = 2
      end if
    else
      if (debug) write(u6,*) ' Field=',line(ifirst:jch-1)
      if (ic == 1) then
        string = line(ifirst:jch-1)
        if (debug) write(u6,*) ' Field read as string :',string
      else if (ic == 2) then
        if (ifirst > jch-1) then
          ierr = 2
          return
        end if
        if (.not. isitanint_cvb(line(ifirst:jch-1))) then
          ierr = 1
          if (debug) write(u6,*) ' Could not read field as integer.'
          return
        end if
        read(line(ifirst:jch-1),*,iostat=istatus) intx
        if (istatus > 0) then
          ierr = 1
          if (debug) write(u6,*) ' Could not read field as integer.'
          return
        end if
        if (debug) write(u6,*) ' Field read as int :',intx
      else if (ic == 3) then
        if (ifirst > jch-1) then
          ierr = 2
          return
        end if
        if (.not. isitareal_cvb(line(ifirst:jch-1))) then
          ierr = 1
          if (debug) write(u6,*) ' Could not read field as real.'
          return
        end if
        read(line(ifirst:jch-1),*,iostat=istatus) realx
        if (istatus > 0) then
          ierr = 1
          if (debug) write(u6,*) ' Could not read field as real.'
          return
        end if
        if (debug) write(u6,*) ' Field read as real :',realx
      end if
    end if
    exit
  end if
end do
if (ich > lenline) then
  write(u6,*) ' Error in input parsing !'
  call abend_cvb()
end if

return

end subroutine gtany_cvb
