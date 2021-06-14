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

#include "unixinfo.fh"
external StrnLn
integer StrnLn
character*35 print

print = ' '
Len = StrnLn(ProgName)
I1 = max(1,35-Len+1)
print(I1:35) = ProgName(1:35-I1+1)
write(6,'(2A)')       ' Program name      :',print
write(6,'(A,I35)')    ' Process ID        :',pid
write(6,'(A,I35)')    ' Parent process ID :',ppid
write(6,'(A,I35)')    ' Seconds           :',sec
write(6,'(A,I35)')    ' Minutes           :',mins
write(6,'(A,I35)')    ' Hours             :',hour
write(6,'(A,I35)')    ' Day of month      :',mday
write(6,'(A,I29,3A)') ' Month             :',mon,' (',Month(mon),')'
write(6,'(A,I35)')    ' Year              :',year
write(6,'(A,I29,3A)') ' Day of week       :',wday,' (',WeekDay(wday),')'
write(6,'(A,I35)')    ' Day of year       :',yday
write(6,'(A,I35)')    ' Daylight saving ? :',isdst
#ifdef _PLEASE_DELETE_ME_
print = ' '
Len = StrnLn(UserName)
I1 = max(1,35-Len+1)
print(I1:35) = UserName(1:35-I1+1)
write(6,'(2A)') ' User name         :',print
print = ' '
Len = StrnLn(RealName)
do i=1,Len
  j = ichar(RealName(i:i))
  if ((j <= 8) .or. (j >= 160) .or. ((j >= 14) .and. (j <= 31))) RealName(i:i) = '?'
end do
I1 = max(1,35-Len+1)
print(I1:35) = RealName(1:35-I1+1)
write(6,'(2A)') ' Name              :',print
print = ' '
Len = StrnLn(HomeDir)
I1 = max(1,35-Len+1)
print(I1:35) = HomeDir(1:35-I1+1)
write(6,'(2A)') ' Home directory    :',print
print = ' '
Len = StrnLn(Shell)
I1 = max(1,35-Len+1)
print(I1:35) = Shell(1:35-I1+1)
write(6,'(2A)') ' Shell             :',print
print = ' '
Len = StrnLn(MolcasDir)
I1 = max(1,35-Len+1)
print(I1:35) = MolcasDir(1:35-I1+1)
write(6,'(2A)') ' Molcas directory  :',print
#endif

return

end subroutine PrtUnixInfo
