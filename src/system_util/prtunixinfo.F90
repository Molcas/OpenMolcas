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

subroutine PrtUnixInfo()

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: I1, length
character(len=35) :: prt
integer(kind=iwp), external :: StrnLn
#include "unixinfo.fh"

prt = ' '
length = StrnLn(ProgName)
I1 = max(1,35-length+1)
prt(I1:35) = ProgName(1:35-I1+1)
write(u6,'(2A)')       ' Program name      :',prt
write(u6,'(A,I35)')    ' Process ID        :',pid
write(u6,'(A,I35)')    ' Parent process ID :',ppid
write(u6,'(A,I35)')    ' Seconds           :',sec
write(u6,'(A,I35)')    ' Minutes           :',mins
write(u6,'(A,I35)')    ' Hours             :',hour
write(u6,'(A,I35)')    ' Day of month      :',mday
write(u6,'(A,I29,3A)') ' Month             :',mon,' (',Month(mon),')'
write(u6,'(A,I35)')    ' Year              :',year
write(u6,'(A,I29,3A)') ' Day of week       :',wday,' (',WeekDay(wday),')'
write(u6,'(A,I35)')    ' Day of year       :',yday
write(u6,'(A,I35)')    ' Daylight saving ? :',isdst

return

end subroutine PrtUnixInfo
