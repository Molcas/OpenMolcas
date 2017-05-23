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
      Subroutine qEnter (String)
************************************************************************
*                                                                      *
*     Initialize an new timer.                                         *
*                                                                      *
*     calling argument:                                                *
*     String : Type character string, input/output                     *
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

#include "SysCtl.fh"
*
      Character*(*) String
      Integer iU,iUinUse
      common /traceiU/iU,iUinUse
#ifdef _DEBUG_
      Character*8   Token
      Integer       Temp(2)
      Equivalence   (Token,Temp)
#endif
*----------------------------------------------------------------------*
*     Initialize the Common / QueCtl / the first time it is referenced *
*----------------------------------------------------------------------*
#ifdef _DEBUG_TRACE_
        if(iUinUse.eq.0) iU=6
        write(iU,'(a,a)') '>> Enter ',string
#endif
#ifdef _DEBUG_
      Temp(1)=0
      Temp(2)=0
      If ( QueCtl(ipStat).ne.ON ) then
         Call IniQue
      End if
*----------------------------------------------------------------------*
*     read default parameters from Common / QueCtl /                   *
*----------------------------------------------------------------------*
      iW=QueCtl(ipSysOut)
      nQTbl=QueCtl(iplTbl)
      nQStk=QueCtl(iplStk)
      If ( QueCtl(ipTrace).eq.ON ) Then
         Write(iW,*) ' <<< Entering qEnter >>>'
      End If
*----------------------------------------------------------------------*
*     convert timer name to standard format                            *
*----------------------------------------------------------------------*
      Temp(1)=0
      Temp(2)=0
      Call StdFmt(String,Token)
      If ( Token.eq.'        ' ) Token='unknown'
      If ( QueCtl(ipTrace).eq.ON ) Then
         Write(iW,'(A,T20,A)') ' Timer name:',Token
      End If
*----------------------------------------------------------------------*
*     search for the table and the stack entry numbers                 *
*----------------------------------------------------------------------*
#ifdef _I8_
      iQTbl=0
      Do 10 i=1,nQTbl
         If ( QueCtl(ipTbl+(i-1)*4+0).eq.Temp(1) ) iQTbl=i
10    Continue
      iQStk=0
      Do 20 i=1,nQStk
         If ( QueCtl(ipStk+(i-1)*3+0).eq.Temp(1) ) iQStk=i
20    Continue
#else
      iQTbl=0
      Do 10 i=1,nQTbl
         If ( QueCtl(ipTbl+(i-1)*4+0).eq.Temp(1) .and.
     &        QueCtl(ipTbl+(i-1)*4+1).eq.Temp(2)       ) iQTbl=i
10    Continue
      iQStk=0
      Do 20 i=1,nQStk
         If ( QueCtl(ipStk+(i-1)*3+0).eq.Temp(1) .and.
     &        QueCtl(ipStk+(i-1)*3+1).eq.Temp(2)       ) iQStk=i
20    Continue
#endif
*----------------------------------------------------------------------*
*     Check the stack entry number                                     *
*----------------------------------------------------------------------*
      If ( iQStk.ne.0 ) then
         Write (6,*) 'QEnter: iQStk.ne.0'
         Write (6,*) 'QEnter: dublicate enter!'
         Write (6,*) 'String=',String
         Write (6,*) 'iQStk=',iQStk
         Call QTrace()
         Call Abend()
      End If
*----------------------------------------------------------------------*
*     Update the stack counter                                         *
*----------------------------------------------------------------------*
      nQStk=nQStk+1
      If ( nQStk.gt.mxQStk ) Then
         Write (6,*) 'QEnter: nQStk.gt.mxQStk'
         Write (6,*) 'nQStk=',nQStk
         Write (6,*) 'mxQStk=',mxQStk
         Call QTrace()
         Call Abend()
      End If
      QueCtl(iplStk)=nQStk
*----------------------------------------------------------------------*
*     Update the stack                                                 *
*----------------------------------------------------------------------*
      iQStk=nQStk
      icpu=iClock()
      QueCtl(ipStk+(iQStk-1)*3+0)=Temp(1)
#ifndef _I8_
      QueCtl(ipStk+(iQStk-1)*3+1)=Temp(2)
#endif
      QueCtl(ipStk+(iQStk-1)*3+2)=icpu
*----------------------------------------------------------------------*
*     Update the table counter                                         *
*----------------------------------------------------------------------*
      If ( iQTbl.eq.0 ) then
         nQTbl=nQTbl+1
         If ( nQTbl.gt.mxQTbl ) Then
            Write (6,*) 'QEnter: nQTbl.gt.mxQTbl'
            Write (6,*) 'nQTbl=',nQTbl
            Write (6,*) 'mxQTbl=',mxQTbl
            Call QTrace()
            Call Abend()
         End If
         QueCtl(iplTbl)=nQTbl
         iQTbl=nQTbl
      End If
*----------------------------------------------------------------------*
*     Update the timer table                                           *
*----------------------------------------------------------------------*
      QueCtl(ipTbl+(iQTbl-1)*4+0)=Temp(1)
#ifndef _I8_
      QueCtl(ipTbl+(iQTbl-1)*4+1)=Temp(2)
#endif
      QueCtl(ipTbl+(iQTbl-1)*4+2)=QueCtl(ipTbl+(iQTbl-1)*4+2)+1
      If ( QueCtl(ipTrace).eq.ON ) Then
         Write(iW,*) ' <<< Exiting qEnter >>>'
         Write(iW,'(A,T20,I3)') ' Table size:',nQTbl
         Write(iW,'(A,T20,I3)') ' Stack size:',nQStk
         Write(iW,'(A,T20,I3)') ' Last entries:',iQTbl
         Write(iW,'(A,T20,I3)') ' Clock time:',icpu
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
