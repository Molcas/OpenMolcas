************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine UnixInfo(SuperModuleName,ModuleName)
      Character*(*) ModuleName, SuperModuleName
      External StrnLn
      Integer StrnLn
#include "unixinfo.fh"
      Logical IfTest
      Data IfTest/.False./

#ifdef _DEBUG_
      IfTest=.True.
#endif

      ProgName=ModuleName
      SuperName=SuperModuleName
      username=' '
      realname=' '
      homedir=' '
      shell=' '
      molcasdir=' '
#ifndef _DEMO_
      Call UnixInfoC(pid, ppid,
     +  sec,mins,hour,mday,mon,year,wday,yday,isdst,
     +  username,realname,homedir,shell,
     +  molcasdir)
#endif

C  Clear ProgName of directory part:
      LenProgName=StrnLn(ProgName)
      Do 100 Ich=LenProgName,1,-1
      If(ProgName(Ich:Ich).Eq.'/')Then
        Islash=Ich
        GoTo 200
      EndIf
100   Continue
      Islash=0
200   Continue
      Do 300 Ich=1,LenProgName
      If(Ich.Le.LenProgName-Islash)Then
        ProgName(Ich:Ich)=ProgName(Ich+Islash:Ich+Islash)
      Else
        ProgName(Ich:Ich)=' '
      EndIf
300   Continue

      mon=mon+1
      year=year+1900
      if(wday.eq.0)wday=7
      yday=yday+1
      WeekDay(1)='Mon'
      WeekDay(2)='Tue'
      WeekDay(3)='Wed'
      WeekDay(4)='Thu'
      WeekDay(5)='Fri'
      WeekDay(6)='Sat'
      WeekDay(7)='Sun'
      Month(1) ='Jan'
      Month(2) ='Feb'
      Month(3) ='Mar'
      Month(4) ='Apr'
      Month(5) ='May'
      Month(6) ='Jun'
      Month(7) ='Jul'
      Month(8) ='Aug'
      Month(9) ='Sep'
      Month(10)='Oct'
      Month(11)='Nov'
      Month(12)='Dec'
      If(IfTest)Call PrtUnixInfo()
      Return
      End
      Subroutine PrtUnixInfo()
#include "unixinfo.fh"
      External StrnLn
      Integer StrnLn
      Character*35 Print

      Print=' '
      Len=StrnLn(ProgName)
      I1=Max(1,35-Len+1)
      Print(I1:35)=ProgName(1:35-I1+1)
      Write(6,'(2A)')' Program name      :',Print
      Write(6,'(A,I35)')' Process ID        :',pid
      Write(6,'(A,I35)')' Parent process ID :',ppid
      Write(6,'(A,I35)')' Seconds           :',sec
      Write(6,'(A,I35)')' Minutes           :',mins
      Write(6,'(A,I35)')' Hours             :',hour
      Write(6,'(A,I35)')' Day of month      :',mday
      Write(6,'(A,I29,3A)')' Month             :',mon,
     +  ' (',Month(mon),')'
      Write(6,'(A,I35)')' Year              :',year
      Write(6,'(A,I29,3A)')' Day of week       :',wday,
     +  ' (',WeekDay(wday),')'
      Write(6,'(A,I35)')' Day of year       :',yday
      Write(6,'(A,I35)')' Daylight saving ? :',isdst
#ifdef _PLEASE_DELETE_ME_
      Print=' '
      Len=StrnLn(UserName)
      I1=Max(1,35-Len+1)
      Print(I1:35)=UserName(1:35-I1+1)
      Write(6,'(2A)')' User name         :',Print
      Print=' '
      Len=StrnLn(RealName)
      Do i=1,Len
         j=iChar(RealName(i:i))
         If ((j.le.8).or.(j.ge.160).or.((j.ge.14).and.(j.le.31)))
     &      RealName(i:i)='?'
      EndDo
      I1=Max(1,35-Len+1)
      Print(I1:35)=RealName(1:35-I1+1)
      Write(6,'(2A)')' Name              :',Print
      Print=' '
      Len=StrnLn(HomeDir)
      I1=Max(1,35-Len+1)
      Print(I1:35)=HomeDir(1:35-I1+1)
      Write(6,'(2A)')' Home directory    :',Print
      Print=' '
      Len=StrnLn(Shell)
      I1=Max(1,35-Len+1)
      Print(I1:35)=Shell(1:35-I1+1)
      Write(6,'(2A)')' Shell             :',Print
      Print=' '
      Len=StrnLn(MolcasDir)
      I1=Max(1,35-Len+1)
      Print(I1:35)=MolcasDir(1:35-I1+1)
      Write(6,'(2A)')' Molcas directory  :',Print
#endif
      Return
      End
