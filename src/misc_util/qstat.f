************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1993, Markus P. Fuelscher                              *
************************************************************************
      Subroutine qStat (String)
************************************************************************
*                                                                      *
*     Print timing information for all sections bracketed by           *
*     calls to qenter and qexit.                                       *
*                                                                      *
*     calling argument:                                                *
*     String: Type character string, input                             *
*             Timer names for which info is printed                    *
*             (string=' ', all info is printed)                        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1993                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "SysCtl.fh"
*
      Character*(*) String
#ifdef _DEBUG_
      Character*8   Token
      Character*8   Fmt1,Fmt2
      Integer       Temp(2)
      Equivalence   (Token,Temp)
      Integer       LocTbl(mxQTbl)
      Integer       LocStk(mxQStk)
      Integer       Substr(2,mxQTbl)
      Character*8   myToken
      Integer       myTemp(2)
      Equivalence   (myToken,myTemp)
*----------------------------------------------------------------------*
      Thr_fraction=0.01D0
      iThr_times=100
      if(iPrintLevel(-1).lt.3) return
*----------------------------------------------------------------------*
*     Initialize the Common / ErrCtl / the first time it is referenced *
*----------------------------------------------------------------------*
      Temp(1)=0
      Temp(2)=0
      If ( QueCtl(ipStat).ne.ON ) then
         Write (6,*) 'QStat: QueCtl(ipStat).ne.ON'
         Call QTrace()
         Call Abend()
      End if
*----------------------------------------------------------------------*
*     read default parameters from Common / QueCtl /                   *
*----------------------------------------------------------------------*
      iW=QueCtl(ipSysOut)
      nQTbl=QueCtl(iplTbl)
      nQStk=QueCtl(iplStk)
      If ( QueCtl(ipTrace).eq.ON ) Then
         Write(iW,*) ' <<< Entering qStat >>>'
      End If
*----------------------------------------------------------------------*
*     check the table counter                                          *
*----------------------------------------------------------------------*
      If ( nQTbl.le.0 ) then
         Write (6,*) 'QStat: Empty table!'
         If ( QueCtl(ipTrace).eq.ON ) Then
            Write(iW,*) ' <<< Exiting qStat >>>'
         End If
         Return
      End If
*----------------------------------------------------------------------*
*     Cut the input string into separate tokens                        *
*----------------------------------------------------------------------*
      lString=LEN(String)
      nToken=0
      iStart=0
      iEnd=0
      Do 10 i=1,lString
        if( lString.eq.1) then
        iStart=1
        iEnd=0
        else
        If ( i.eq.1 .and. String(i:i).ne.' ' ) then
           iStart=1
        Else if ( String(i:i).ne.' ' .and. String(i-1:i-1).eq.' ' )then
           iStart=i
        Else if ( String(i:i).eq.' ' .and. String(i-1:i-1).ne.' ' )then
           iEnd=i-1
        End If
        endif
        If (i.eq.lString.and.String(lString:lString).ne.' ')iEnd=lString
        If ( iStart.ne.0 .and. iEnd.ne.0 ) then
          nToken=nToken+1
          SubStr(1,nToken)=iStart
          SubStr(2,nToken)=iEnd
          iStart=0
          iEnd=0
        End If
10    Continue
*----------------------------------------------------------------------*
*     Make o local copy of the timer table and stack                   *
*----------------------------------------------------------------------*
      Do 20 i=1,nQTbl
         LocTbl(i)=QueCtl(ipTbl+(i-1)*4+3)
20    Continue
      Do 25 i=1,nQStk
         LocStk(i)=QueCtl(ipStk+(i-1)*3+2)
25    Continue
*----------------------------------------------------------------------*
*     Update the local stack using the current time                    *
*----------------------------------------------------------------------*
      icpu=iClock()
      Do 30 i=nQStk,1,-1
         LocStk(i)=icpu-LocStk(i)
30    Continue
      Do 32 i=1,nQStk-1
         LocStk(i)=LocStk(i)-LocStk(i+1)
32    Continue
*----------------------------------------------------------------------*
*     Update the local timer table using the current time              *
*----------------------------------------------------------------------*
      icpu=iClock()
      Do 35 i=nQStk,1,-1
         Do 37 j=1,nQTbl
            If (QueCtl(ipStk+(i-1)*3+0).eq.QueCtl(ipTbl+(j-1)*4+0) .and.
     &          QueCtl(ipStk+(i-1)*3+1).eq.QueCtl(ipTbl+(j-1)*4+1) )
     &          LocTbl(j)=LocTbl(j)+LocStk(i)
37       Continue
35    Continue
*----------------------------------------------------------------------*
*     Finally, print the timings                                       *
*----------------------------------------------------------------------*
      Fmt1='(2X,A)'
      Fmt2='(2X,'
      icput=0
      Do 40 i=1,nQTbl
         icput=icput+LocTbl(i)
40    Continue
      icput=Max(1,icput)
      Write(iW,*)
      Call CollapseOutput(1,'Statistics and timing')
      Write(iW,'(3X,A)')    '---------------------'
      Write(iW,*)
      If ( nToken.eq.0 ) then
         Write(iW,Fmt1)
     &   'Name     Num. of    tot cpu time     fraction   '
         Write(iW,Fmt1)
     &   '          calls        [sec.]                   '
         Write(iW,Fmt1)
     &   '- - - - - - - - - - - - - - - - - - - - - - - - '
         Do 50 i=1,nQTbl
            ttot=DBLE(LocTbl(i))/DBLE(ClkInc)
            fraction=LocTbl(i)/DBLE(icput)
            If (fraction.gt.Thr_fraction .or.
     &          QueCtl(ipTbl+(i-1)*4+2).gt.iThr_times) Then
#ifdef _I8_
            myTemp(1)=QueCtl(ipTbl+(i-1)*4+0)

            Write(iW,Fmt2//' A8,I6,2F15.2)')
c     &           QueCtl(ipTbl+(i-1)*4+0),
     &            myToken,
     &           QueCtl(ipTbl+(i-1)*4+2),
     &           ttot,fraction
#else
            myTemp(1)=QueCtl(ipTbl+(i-1)*4+0)
            myTemp(2)=QueCtl(ipTbl+(i-1)*4+1)
            Write(iW,Fmt2//' A8,I6,2F15.2)')
c     &           QueCtl(ipTbl+(i-1)*4+0),
c     &           QueCtl(ipTbl+(i-1)*4+1),
     &           myToken,
     &           QueCtl(ipTbl+(i-1)*4+2),
     &           ttot,fraction
#endif
            End If
50       Continue
      Else
         Write(iW,Fmt1)
     &   'Name     Num. of    tot cpu time     fraction   '
         Write(iW,Fmt1)
     &   '          calls        [sec.]                   '
         Write(iW,Fmt1)
     &   '- - - - - - - - - - - - - - - - - - - - - - - - '
         Do 60 i=1,nToken
            Temp(1)=0
            Temp(2)=0
            Call StdFmt(String(SubStr(1,i):SubStr(2,i)),Token)
            Do 65 j=1,nQtbl
#ifdef _I8_
               If ( QueCtl(ipTbl+(j-1)*4+0).eq.Temp(1) ) then
#else
               If ( QueCtl(ipTbl+(j-1)*4+0).eq.Temp(1) .and.
     &              QueCtl(ipTbl+(j-1)*4+1).eq.Temp(2)       ) then
#endif
                    ttot=DBLE(LocTbl(j))/DBLE(ClkInc)
                    fraction=LocTbl(j)/DBLE(icput)
               If (fraction.gt.Thr_fraction .or.
     &             QueCtl(ipTbl+(i-1)*4+2).gt.iThr_times) Then
#ifdef _I8_
            myTemp(1)=QueCtl(ipTbl+(j-1)*4+0)

                    Write(iW,Fmt2//' A8,I6,2F15.2)')
c     &                    QueCtl(ipTbl+(j-1)*4+0),
     &                    myToken,
     &                    QueCtl(ipTbl+(j-1)*4+2),
     &                    ttot,fraction
#else
            myTemp(1)=QueCtl(ipTbl+(j-1)*4+0)
            myTemp(2)=QueCtl(ipTbl+(j-1)*4+1)
                    Write(iW,Fmt2//'A8,I6,2F15.2)')
c     &                    QueCtl(ipTbl+(j-1)*4+0),
c     &                    QueCtl(ipTbl+(j-1)*4+1),
     &                     myToken,
     &                    QueCtl(ipTbl+(j-1)*4+2),
     &                    ttot,fraction
#endif
               End If
               End If
65          Continue
60       Continue
      End If
      ttot=DBLE(icput)/DBLE(ClkInc)
      fraction=1.0
      Write(iW,Fmt1) '- - - - - - - - - - - - - - - - - - - - - - - - '
      Write(iW,Fmt2//'A,6X,2F15.2)') 'Total   ',ttot,fraction
      Write(iW,Fmt1) '- - - - - - - - - - - - - - - - - - - - - - - - '
      Call CollapseOutput(0,'Statistics and timing')
      Write(iW,*)
      If ( QueCtl(ipTrace).eq.ON ) Then
         Write(iW,*) ' <<< Exiting qStat >>>'
      End If
*----------------------------------------------------------------------*
*     exit                                                             *
*----------------------------------------------------------------------*
#else
c Avoid unused argument warnings
      If (.False.) Call Unused_character(String)
#endif
      Return
      End
