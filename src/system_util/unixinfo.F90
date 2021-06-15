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

subroutine UnixInfo(SuperModuleName,ModuleName)

use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: SuperModuleName, ModuleName
integer(kind=iwp) :: Ich, Islash, LenProgName
integer(kind=iwp), external :: StrnLn
logical(kind=iwp) :: IfTest = .false.
#include "unixinfo.fh"

#ifdef _DEBUGPRINT_
IfTest = .true.
#endif

ProgName = ModuleName
SuperName = SuperModuleName
username = ' '
realname = ' '
homedir = ' '
shell = ' '
molcasdir = ' '
#ifndef _DEMO_
call UnixInfoC(pid,ppid,sec,mins,hour,mday,mon,year,wday,yday,isdst,username,realname,homedir,shell,molcasdir)
#endif

!  Clear ProgName of directory part:
LenProgName = StrnLn(ProgName)
do Ich=LenProgName,1,-1
  if (ProgName(Ich:Ich) == '/') then
    Islash = Ich
    goto 200
  end if
end do
Islash = 0
200 continue
do Ich=1,LenProgName
  if (Ich <= LenProgName-Islash) then
    ProgName(Ich:Ich) = ProgName(Ich+Islash:Ich+Islash)
  else
    ProgName(Ich:Ich) = ' '
  end if
end do

mon = mon+1
year = year+1900
if (wday == 0) wday = 7
yday = yday+1
WeekDay(1) = 'Mon'
WeekDay(2) = 'Tue'
WeekDay(3) = 'Wed'
WeekDay(4) = 'Thu'
WeekDay(5) = 'Fri'
WeekDay(6) = 'Sat'
WeekDay(7) = 'Sun'
Month(1) = 'Jan'
Month(2) = 'Feb'
Month(3) = 'Mar'
Month(4) = 'Apr'
Month(5) = 'May'
Month(6) = 'Jun'
Month(7) = 'Jul'
Month(8) = 'Aug'
Month(9) = 'Sep'
Month(10) = 'Oct'
Month(11) = 'Nov'
Month(12) = 'Dec'
if (IfTest) call PrtUnixInfo()

return

end subroutine UnixInfo
