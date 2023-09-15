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
logical debug
logical isitanint_cvb, isitareal_cvb
external isitanint_cvb, isitareal_cvb
#include "luinp_cvb.fh"
#include "rdline.fh"
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
1100 jch = jch+1
    if ((ilv(jch) == 0) .and. (jch /= lenline+1)) goto 1100
    ifirst = ich+1
    if (ich == 1) ifirst = 1
    ! Special character strings to signify empty field?
    do iempty=1,nempty
      if (line(ifirst:jch-1) == empty(iempty)(1:len_trim_cvb(empty(iempty)))) goto 1130
    end do
    if (debug) write(6,*) ' Field=',line(ifirst:jch-1)
    if (ic == 1) then
      string = line(ifirst:jch-1)
      if (debug) write(6,*) ' Field read as string :',string
    else if (ic == 2) then
      if (ifirst > jch-1) then
        ierr = 2
        return
      end if
      if (.not. isitanint_cvb(line(ifirst:jch-1))) goto 1150
      read(line(ifirst:jch-1),*,err=1150) int
      if (debug) write(6,*) ' Field read as int :',int
    else if (ic == 3) then
      if (ifirst > jch-1) then
        ierr = 2
        return
      end if
      if (.not. isitareal_cvb(line(ifirst:jch-1))) goto 1150
      read(line(ifirst:jch-1),*,err=1150) real
      if (debug) write(6,*) ' Field read as real :',real
    end if
    return
1130 continue
    ! "Empty" field:
    if (ic == 1) then
      string = ' '
    else
      ierr = 2
    end if
    return
  end if
end do
write(6,*) ' Error in input parsing !'
call abend_cvb()
1150 ierr = 1
if (debug) then
  if (ic == 2) then
    write(6,*) ' Could not read field as integer.'
  else
    write(6,*) ' Could not read field as real.'
  end if
end if

return

end subroutine gtany_cvb
