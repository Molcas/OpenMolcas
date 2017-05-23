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
      Subroutine qExit (String)
************************************************************************
*                                                                      *
*     Close a timer.                                                   *
*                                                                      *
*     calling argument:                                                *
*     String : Type character string, input/output                     *
*              Timer name.                                             *
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
#ifdef _DEBUG_
      Character*8   Token
      Integer       Temp(2)
      Equivalence   (Token,Temp)
*----------------------------------------------------------------------*
*     Initialize the Common / ErrCtl / the first time it is referenced *
*----------------------------------------------------------------------*
      Temp(1)=0
      Temp(2)=0
      If ( QueCtl(ipStat).ne.ON ) then
         Write (6,*) 'QExit: QueCtl(ipStat).ne.ON'
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
         Write(iW,*) ' <<< Entering qExit >>>'
         Write(iW,'(A,T20,I3)') ' Table size:',nQTbl
         Write(iW,'(A,T20,I3)') ' Stack size:',nQStk
      End If
*----------------------------------------------------------------------*
*     check the stack counter                                          *
*----------------------------------------------------------------------*
      If ( nQStk.le.0 ) Then
         Write (6,*) 'QExit: String=',String
         Write (6,*) 'QExit: nQStk.le.0'
         Call QTrace()
         Call Abend()
      End if

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
      If ( iQStk.eq.0 ) then
         Write (6,*) 'QExit: iQStk.eq.0'
         Write (6,*) 'String=',String
         Write (6,*) 'Subroutine not entered!'
         Call QTrace()
         Call Abend()
      End If
      If ( iQStk.ne.nQStk ) then
         Write (6,*) 'QExit: iQStk.ne.nQStk'
         Write (6,*) 'String=',String
         Write (6,*) 'Exit out of sequence!'
         Call QTrace()
         Call Abend()
      End If
*----------------------------------------------------------------------*
*     Check the table entry number                                     *
*----------------------------------------------------------------------*
      If ( iQTbl.eq.0 ) then
         Write (6,*) 'QExit: iQTbl.eq.0'
         Write (6,*) 'String=',String
         Call QTrace()
         Call Abend()
      End If
*----------------------------------------------------------------------*
*     Update the stack                                                 *
*----------------------------------------------------------------------*
      icpu=iClock()
      jcpu=icpu-QueCtl(ipStk+(iQStk-1)*3+2)
      Do 30 i=nQStk-1,1,-1
         QueCtl(ipStk+(i-1)*3+2)=QueCtl(ipStk+(i-1)*3+2)+jcpu
30    Continue
*----------------------------------------------------------------------*
*     Update the timer table                                           *
*----------------------------------------------------------------------*
      QueCtl(ipTbl+(iQTbl-1)*4+3)=QueCtl(ipTbl+(iQTbl-1)*4+3)+jcpu
*----------------------------------------------------------------------*
*     Update the stack counter                                         *
*----------------------------------------------------------------------*
      nQStk=nQStk-1
      QueCtl(iplStk)=nQStk
      If ( QueCtl(ipTrace).eq.ON ) Then
         Write(iW,*) ' <<< Exiting qExit >>>'
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
