!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine getSeed(iseed)

use UnixInfo, only: ProgName
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: iseed
integer(kind=iwp) :: days, hours, i, minutes, seconds
character(len=72) :: Line

! Externally defined seed
call getenvf('MOLCAS_RANDOM_SEED',Line)
if (Line /= ' ') then
  read(Line,*) iseed
  return
end if

! Somewhat reproducible if inside verification
call getenvf('MOLCAS_TEST',Line)
if (Line /= ' ') then
  call getenvf('MOLCAS_ITER',Line)
  read(Line,*) iseed
  call getenvf('MOLCAS_PRINT',Line)
  do i=1,len_trim(Line)
    iseed = iseed+ichar(Line(i:i))
  end do
  Line = trim(ProgName)
  do i=1,len_trim(Line)
    iseed = iseed+ichar(Line(i:i))
  end do
  return
end if

! Default: based on time and project name
call datimx(Line)
read(Line,'(8x,i2,1x,i2,1x,i2,1x,i2)') days,hours,minutes,seconds
iseed = int(((days*24+hours)*60+minutes)*60+seconds,kind(iseed))
call getenvf('Project',Line)
do i=1,len_trim(Line)
  iseed = iseed+ichar(Line(i:i))
end do

return

end subroutine getSeed
