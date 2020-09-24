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
* Copyright (C) 1993, Johan Lorentzon                                  *
*               1993, Jeppe Olsen                                      *
*               1993, Markus P. Fuelscher                              *
************************************************************************
      Subroutine Rd2Int(iPL)
************************************************************************
*                                                                      *
*     Read header of the two-electron integral file                    *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     J. Lorentzon, J. Olsen and M.P. Fuelscher                        *
*     University of Lund, Sweden, 1993                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)

#include "Input.fh"
      Integer nSymX,nBasX(mxSym)
      Logical SqSym
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
      iRc=-1
      Call GetOrd(iRc,SqSym,nSymX,nBasX,nSkip)
      If ( iRc.ne.0 ) Then
         Write (6,*) 'Rd2Int: Error reading ORDINT'
         Call QTrace
         Call Abend()
      End If
      If (iPL.ge.2) Then
      If(SqSym)write(6,*)'OrdInt status: squared'
      If(.not.SqSym)write(6,*)'OrdInt status: non-squared'
      End If
      If ( nSymX.ne.nSym ) Then
         Write (6,*) 'Rd2Int: nSymX.ne.nSym'
         Write (6,*) 'nSymX,nSym=',nSymX,nSym
         Call QTrace
         Call Abend()
      End If
      Do 10 iSym=1,nSym
         If ( nBas(iSym).ne.nBasX(iSym) ) Then
            Write (6,*) 'Rd2Int: nBas(iSym).ne.nBasX(iSym)'
            Write (6,*) 'nBas(iSym),nBasX(iSym)=',
     &                   nBas(iSym),nBasX(iSym)
            Call QTrace
            Call Abend()
         End If
10    Continue
      ntSkip=0
      Do 20 iSym=1,nSym
         ntSkip=ntSkip+nSkip(iSym)
20    Continue
      If ( ntSkip.ne.0 ) Then
         Write (6,*) 'Rd2Int: ntSkip.ne.0'
         Write (6,*) 'ntSkip=',ntSkip
         Call QTrace
         Call Abend()
      End If
      If (.not.SqSym.and..not.TimeDep) Then
         CASINT=.true.
      Else
         CASINT=.False.
      End If
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
      Return
      End
