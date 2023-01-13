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

module UnixInfo

use Definitions, only: iwp

implicit none
private

integer(kind=iwp), protected :: hour, isdst, mday, mins, mon, pid, ppid, sec, wday, yday, year
character(len=256), protected :: HomeDir, MolcasDir, ProgName, RealName, Shell, SuperName, UserName
character(len=*), parameter :: Month(12) = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'], &
                               Weekday(7) = ['Mon','Tue','Wed','Thu','Fri','Sat','Sun']

public :: HomeDir, hour, init_UnixInfo, isdst, mday, mins, MolcasDir, mon, Month, pid, ppid, ProgName, RealName, sec, Shell, &
          SuperName, UserName, wday, Weekday, yday, year

contains

subroutine init_UnixInfo(SuperModuleName,ModuleName)

  character(len=*), intent(in) :: SuperModuleName, ModuleName
  integer(kind=iwp) :: Ich, Islash, LenProgName
  integer(kind=iwp), external :: StrnLn
  logical(kind=iwp) :: IfTest = .false.
  interface
    subroutine UnixInfoC(pid,ppid,sec,min_,hour,mday,mon,year,wday,yday,isdst,username,realname,homedir,shell,molcasdir) &
               bind(C,name='unixinfoc_')
      use, intrinsic :: iso_c_binding, only: c_char
      use Definitions, only: MOLCAS_C_INT
      integer(kind=MOLCAS_C_INT) :: pid, ppid, sec, min_, hour, mday, mon, year, wday, yday, isdst
      character(kind=c_char) :: username(*), realname(*), homedir(*), shell(*), molcasdir(*)
    end subroutine UnixInfoC
  end interface

# ifdef _DEBUGPRINT_
  IfTest = .true.
# endif

  ProgName = ModuleName
  SuperName = SuperModuleName
  UserName = ' '
  RealName = ' '
  HomeDir = ' '
  Shell = ' '
  MolcasDir = ' '
# ifndef _DEMO_
  call UnixInfoC(pid,ppid,sec,mins,hour,mday,mon,year,wday,yday,isdst,UserName,RealName,HomeDir,Shell,MolcasDir)
# endif

  ! Clear ProgName of directory part:
  LenProgName = StrnLn(ProgName)
  Islash = 0
  do Ich=LenProgName,1,-1
    if (ProgName(Ich:Ich) == '/') then
      Islash = Ich
      exit
    end if
  end do
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
  if (IfTest) call PrtUnixInfo()

  return

end subroutine init_UnixInfo

subroutine PrtUnixInfo()

  use Definitions, only: u6

  character(len=35) :: prt

  prt = adjustr(ProgName(1:35))
  write(u6,'(2A)')       ' Program name      :',prt
  write(u6,'(A,I35)')    ' Process ID        :',pid
  write(u6,'(A,I35)')    ' Parent process ID :',ppid
  write(u6,'(A,I35)')    ' Seconds           :',sec
  write(u6,'(A,I35)')    ' Minutes           :',mins
  write(u6,'(A,I35)')    ' Hours             :',hour
  write(u6,'(A,I35)')    ' Day of month      :',mday
  write(u6,'(A,I29,3A)') ' Month             :',mon,' (',Month(mon),')'
  write(u6,'(A,I35)')    ' Year              :',year
  write(u6,'(A,I29,3A)') ' Day of week       :',wday,' (',Weekday(wday),')'
  write(u6,'(A,I35)')    ' Day of year       :',yday
  write(u6,'(A,I35)')    ' Daylight saving ? :',isdst

  return

  end subroutine PrtUnixInfo

end module UnixInfo
