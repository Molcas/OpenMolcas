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
* Copyright (C) Roland Lindh                                           *
*               Giovanni Ghigo                                         *
************************************************************************
      SubRoutine Rd_UDIC(nLines,iInt,nFix,nRowH)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
*             University of Lund, SWEDEN                               *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Character*120 Temp
      character*16 filnam
*
      iRout = 27
      iPrint = nPrint(iRout)
*
      Lu=6
*
      Lu_UDIC=91
      filnam='UDIC'
      call molcas_open(Lu_UDIC,filnam)
c      Open(Lu_UDIC,File=filnam,Form='FORMATTED',Status='OLD')
      Rewind(Lu_UDIC)
*
*     Find begining of definitions of internal coordinates
*
      Do iLines = 1, nLines
         Read(Lu_UDIC,'(A)') Temp
         Call UpCase(Temp)
         If (Temp(1:4).eq.'VARY') Go To 100
      End Do
      Call WarningMessage(2,' No internal coordinates are defined!')
      Call Quit_OnUserError()
*
 100  Continue
      iInt = 0
      nFix = 0
      nRowH = 0 ! Number of Rows of Hessian Numerically estimated
      Do jLines = iLines+1, nLines
         Read(Lu_UDIC,'(A)') Temp
         Call UpCase(Temp)
         If (Temp(1:3).eq.'FIX') Go To 200
         If (Temp(1:4).eq.'ROWH') then
           kLines = jLines
           Go To 300
         EndIf
*------- Do not count line if continuation character
         If (Index(Temp,'&').eq.0) iInt=iInt+1
      End Do
      Go To 400
*
 200  Continue
      Do kLines = jLines+1, nLines
         Read(Lu_UDIC,'(A)') Temp
         Call UpCase(Temp)
         If (Temp(1:4).eq.'ROWH') Go To 300
*------- Do not count line if continuation character
         If (Index(Temp,'&').eq.0) nFix=nFix+1
      End Do
 300  Do lLines = kLines+1, nLines
         Read(Lu_UDIC,'(A)') Temp
         Call UpCase(Temp)
*------- Do not count line if continuation character
         If (Index(Temp,'&').eq.0) nRowH=nRowH+1
      End Do
 400  Continue
*
      Close(Lu_UDIC)
      Return
      End


      SubRoutine Rd_UDIC_RowH(nInter,nRowH,mRowH)
************************************************************************
*                                                                      *
* Object: Reading the Internal Coordinates required for Numerical      *
*         estimation of single rows and columns of Hessian             *
* Called from: PrePro when nRowH.GT.0                                  *
* Author: Giovanni Ghigo, University of Torino, Italy                  *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Character*8 Labels(nInter)
      Character*8 cLbl
      Character*120 Temp
      Character*16 filnam
      Integer mRowH(10)
*
*
* nMaxRowH is the maximum number of internal coordinates that can be
* specified for the numerical estimation of rows (and column) of the
* hessian matrix. The mRowH(10) is defined everywhere is used. There
* are no limitation, in principle, to nMaxRowH.
*
      nMaxRowH = 10
*
      Lu=6
      Lu_UDIC=91
      filnam='UDIC'
      call molcas_open(Lu_UDIC,filnam)
      Rewind(Lu_UDIC)
      Do iRowH = 1, nMaxRowH
         mRowH(iRowH) = 0
      EndDo
      If (nRowH.GT.nMaxRowH) then
         Call WarningMessage(2,'Error in rd_udic')
         Write (Lu,*) '***************************************'
         Write (Lu,*) ' ERROR: Too many coordinates in ROWH ! '
         Write (Lu,*) '        Max ',nMaxRowH
         Write (Lu,*) '***************************************'
         Call Quit_OnUserError()
      EndIf
*
*     Find begining of definitions of internal coordinates
*
 10   Read(Lu_UDIC,'(A)') Temp
      Call UpCase(Temp)
      If (Temp(1:4).ne.'VARY') Go To 10
      Do iLines = 1, nInter
 20      Read(Lu_UDIC,'(A)') Temp
         Call UpCase(Temp)
         If (Temp(1:3).eq.'FIX' ) Go To 20
         cLbl = '        '
         Do j = 1, 120
            If (Temp(j:j).EQ.' ') GoTo 30
            cLbl(j:j) = Temp(j:j)
         Enddo
 30      Labels(iLines) = cLbl
 35      If (Index(Temp,'&').eq.0) GoTo 37
         Read(Lu_UDIC,'(A)') Temp
         GoTo 35
 37   Continue
      EndDo
*
      Read(Lu_UDIC,'(A)') Temp ! Skip ROWH
      Do iRowH = 1 , nRowH
         Read(Lu_UDIC,'(A)') Temp
         Call UpCase(Temp)
         cLbl(1:8) = Temp(1:8)
         Do kLines = 1, nInter
            If (cLbl.eq.Labels(kLines)) then
               mRowH(iRowH) = kLines
               GoTo 40
            EndIf
         EndDo
         Call WarningMessage(2,'Error in rd_udic')
         Write (Lu,*) '**********************************************'
         Write (Lu,*) ' ERROR: Undefined internal ROWH coordinate in '
         Write (Lu,*) ' ',Temp(1:60)
         Write (Lu,*) '**********************************************'
         Call Quit_OnUserError()
 40      Continue
      End Do
      Close(Lu_UDIC)
      Return
      End
