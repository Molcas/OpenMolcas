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

subroutine gtany_cvb(string,int,real,ic,ifield,ierr)

implicit real*8(a-h,o-z)
logical debug, done
logical isitanint_cvb, isitareal_cvb
external isitanint_cvb, isitareal_cvb
#include "luinp_cvb.fh"
#include "rdline.fh"
integer istatus
parameter(nblank=2,nempty=1)
character*(*) string
character*2 empty(nempty)
save empty
data empty/'--'/
data debug/.false./

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
      if (line(ifirst:jch-1) == empty(iempty)(1:len_trim_cvb(empty(iempty)))) then
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
      if (debug) write(6,*) ' Field=',line(ifirst:jch-1)
      if (ic == 1) then
        string = line(ifirst:jch-1)
        if (debug) write(6,*) ' Field read as string :',string
      else if (ic == 2) then
        if (ifirst > jch-1) then
          ierr = 2
          return
        end if
        if (.not. isitanint_cvb(line(ifirst:jch-1))) then
          ierr = 1
          if (debug) write(6,*) ' Could not read field as integer.'
          return
        end if
        read(line(ifirst:jch-1),*,iostat=istatus) int
        if (istatus > 0) then
          ierr = 1
          if (debug) write(6,*) ' Could not read field as integer.'
          return
        end if
        if (debug) write(6,*) ' Field read as int :',int
      else if (ic == 3) then
        if (ifirst > jch-1) then
          ierr = 2
          return
        end if
        if (.not. isitareal_cvb(line(ifirst:jch-1))) then
          ierr = 1
          if (debug) write(6,*) ' Could not read field as real.'
          return
        end if
        read(line(ifirst:jch-1),*,iostat=istatus) real
        if (istatus > 0) then
          ierr = 1
          if (debug) write(6,*) ' Could not read field as real.'
          return
        end if
        if (debug) write(6,*) ' Field read as real :',real
      end if
    end if
    exit
  end if
end do
if (ich > lenline) then
  write(6,*) ' Error in input parsing !'
  call abend_cvb()
end if

return

end subroutine gtany_cvb
