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
* Copyright (C) 2001, Roland Lindh                                     *
*               Sergey Gusarov                                         *
************************************************************************
      SubRoutine GetPAM(lUnit,ipExp,ipCff,nExp,nBasis,MxShll,iShll,
     &                  BLine,
     &                  ipPAM2xp,ipPAM2cf,nPAM2,DInf,nDInf)
************************************************************************
*                                                                      *
*    Objective: To read potential information, for DMFT calculations   *
*               This means that we read (from input stream)            *
*               the potential terms (PAM)                              *
*                                                                      *
* Called from: GetBs                                                   *
*                                                                      *
* Calling    : RecPrt, Rdbsl                                           *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*                                                                      *
*     Modified: Sergey Gusarov SPb State Ubiv, Russia.                 *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "print.fh"
      Real*8 DInf(nDInf)
#include "real.fh"
      Character*80 BLine
      Character*180 Line, Get_Ln
*     External Get_Ln
      Integer ipExp(MxShll), ipCff(MxShll), nExp(MxShll), nBasis(MxShll)
      Logical test
      data test /.False./
c     data test /.True./
      iRout=6
      iPrint = nPrint(iRout)
      Call QEnter('GetPAM')
*
      if (test) Write (6,*) ' Reading PAM potencials'
      iStrt = ipExp(iShll+1)
      ipPAM2xp = iStrt
*     Read(Line,*) nPAM2
      Line=Get_Ln(lUnit)
      If (Index(Line,'PAM').eq.0) Then
         Call WarningMessage(2,
     &              'ERROR: Keyword PAM expected, offending line : '
     &            //Line)
         Call Quit_OnUserError()
      Endif
      Line=Get_Ln(lUnit)
      Call Get_i(1,nPAM2,1)
      Do iPAM_Ang=0, nPAM2
         Line=Get_Ln(lUnit)
         Call Get_i(1,nPAM2Prim,1)
         Call Get_i(2,nPAM2Bas,1)
         DInf(iStrt)=nPAM2Prim
         DInf(iStrt+1)=nPAM2Bas
         iStrt=iStrt+2
         iEnd = iStrt + nPAM2Prim - 1
*
*   Read exponents
*
         If (nPAM2Prim.gt.0) then
            Call read_v(lUnit,DInf,iStrt,iEnd,1,Ierr)
            If (Ierr.ne.0) Then
               Call WarningMessage(2,
     &                     'GetBS: Error reading GPA exponents')
               Call Abend()
            End If
         End If
*
*   Read coefficents
*
      iStrt = iEnd + 1
         if (iPAM_Ang.eq.0) ipPAM2cf = iStrt
         iEnd  = iStrt + nPAM2Prim*nPAM2Bas - 1
         Do 288 iPrim = 0, nPAM2Prim - 1
            Call Read_v(lUnit,DInf,iStrt+iPrim,iEnd,nPAM2Prim,ierr)
            If(ierr.ne.0) Then
               Call WarningMessage(2,
     &                         'GetBS: Error in reading GPA!!!')
               Call Abend()
            End if
 288     Continue
         iStrt = iEnd + 1
      End Do
      ipExp(iShll+1) = iEnd + 1
*
      Call QExit('GetPAM')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer_array(ipCff)
         Call Unused_integer_array(nExp)
         Call Unused_integer_array(nBasis)
         Call Unused_character(BLine)
      End If
      End
*
