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
      Subroutine qTrace
************************************************************************
*                                                                      *
*     Print trace back of calling routines                             *
*                                                                      *
*     calling arguments: none                                          *
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
#ifdef _DEBUGPRINT_
      Character*8   Proc1,Proc2
      Integer       iProc1(2),iProc2(2)
      Equivalence   (Proc1,iProc1),(Proc2,iProc2)
*----------------------------------------------------------------------*
*     Initialize the Common / ErrCtl / the first time it is referenced *
*----------------------------------------------------------------------*
      If ( QueCtl(ipStat).ne.ON ) then
         Write (6,*) 'QTrace: QueCtl(ipStat).ne.ON'
         Call Abend()
      End if
*----------------------------------------------------------------------*
*     read default parameters from Common / QueCtl /                   *
*----------------------------------------------------------------------*
      iW=QueCtl(ipSysOut)
      nQTbl=QueCtl(iplTbl)
      nQStk=QueCtl(iplStk)
      If ( QueCtl(ipTrace).eq.ON ) Then
         Write(iW,*) ' <<< Entering qTrace >>>'
      End If
*----------------------------------------------------------------------*
*     check the stack counter                                          *
*----------------------------------------------------------------------*
      If ( nQStk.le.0 ) Then
         Write (6,*) 'QTrace: no routines to trace!'
         Return
      End If
*----------------------------------------------------------------------*
*     dump the current content remark system                           *
*----------------------------------------------------------------------*
c      Call RemDmp(iW,'FIFO')
*----------------------------------------------------------------------*
*     dump the current content of the stack                            *
*----------------------------------------------------------------------*
      If ( nQStk.ge.1 ) then
        Write(iW,*)
        Write(iW,*) ' Calling history'
        Write(iW,*) ' ---------------'
        Write(iW,*)
      End If
      iProc1(1)=QueCtl(ipStk+(nQStk-1)*3+0)
      iProc1(2)=QueCtl(ipStk+(nQStk-1)*3+1)
      Write(iW,'(2X,2A)') 'last entry: ',Proc1
      Do i=nQStk,2,-1
        iProc1(1)=QueCtl(ipStk+(i-1)*3+0)
        iProc1(2)=QueCtl(ipStk+(i-1)*3+1)
        iProc2(1)=QueCtl(ipStk+(i-2)*3+0)
        iProc2(2)=QueCtl(ipStk+(i-2)*3+1)
        Write(iW,'(2X,3A)') Proc1,' called by ',Proc2
      End Do
      iProc1(1)=QueCtl(ipStk+0)
      iProc1(2)=QueCtl(ipStk+1)
      Write(iW,'(2X,2A)') 'root:       ',Proc1
      Write(iW,*)
      If ( QueCtl(ipTrace).eq.ON ) Then
         Write(iW,*) ' <<< Exiting qTrace >>>'
      End If
*----------------------------------------------------------------------*
*     exit                                                             *
*----------------------------------------------------------------------*
#endif
      Return
      End
