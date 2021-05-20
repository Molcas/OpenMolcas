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
* Copyright (C) 1991, Roland Lindh                                     *
************************************************************************
      SubRoutine NxtWrd(Line,iF,iE)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
* Called from: DefInt                                                  *
*                                                                      *
* Calling    : None                                                    *
*                                                                      *
*     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
*             University of Lund, SWEDEN                               *
*             May '91                                                  *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Character*(*) Line
*
      nChar = Len(Line)
*     Find first non-blank character
 10   Continue
      If (iF.eq.0.or.iF.gt.nChar) Then
         Call WarningMessage(2,'NxtWrd: iF.eq.0.or.iF.gt.nChar')
         Write (6,*) 'nChar=',nChar
         Write (6,*) 'iF,iE=',iF,iE
         Call Abend()
      End If
      If (Line(iF:iF).eq.' ') Then
         iF = iF + 1
         If (iF.ge.nChar) Then
             iF = nChar
            iE=-1
            Return
         End If
         Go To 10
      End If
*     Find the end of the present word
      iE = iF + 1
 20   Continue
      If (Line(iE:iE).ne.' ') Then
         iE = iE + 1
         If (iE.gt.nChar) Then
            iE=nChar
            Return
         End If
         Go To 20
      Else
         iE = iE - 1
      End If
*
      Return
      End
