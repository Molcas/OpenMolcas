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
* Copyright (C) 2001,2020,  Roland Lindh                               *
*               Sergey Gusarov                                         *
************************************************************************
      SubRoutine GetPAM(lUnit,iCnttp)
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
      Use Basis_Info
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "real.fh"
#include "stdalloc.fh"
      Real*8, Allocatable:: Array(:)
      Character*180 Line, Get_Ln
*     External Get_Ln
      Integer nPAM2
      Logical test
#ifdef _DEBUG_
      data test /.True./
#else
      data test /.False./
#endif
*
      if (test) Write (6,*) ' Reading PAM potencials'
      nArray=10000
      Call mma_Allocate(Array,nArray,Label='Array')
*
      iStrt = 1
      Line=Get_Ln(lUnit)
      If (Index(Line,'PAM').eq.0) Then
         Call WarningMessage(2,
     &              'ERROR: Keyword PAM expected, offending line : '
     &            //Line)
         Call Quit_OnUserError()
      Endif
      Line=Get_Ln(lUnit)
      Call Get_i1(1,nPAM2)
      dbsc(iCnttp)%nPAM2=nPAM2
      Do iPAM_Ang=0, nPAM2
         Line=Get_Ln(lUnit)
         Call Get_i1(1,nPAM2Prim)
         Call Get_i1(2,nPAM2Bas)
         Array(iStrt)=DBLE(nPAM2Prim)
         Array(iStrt+1)=DBLE(nPAM2Bas)
         iStrt=iStrt+2
         iEnd = iStrt + nPAM2Prim - 1
*
*   Read exponents
*
         If (nPAM2Prim.gt.0) then
            Call read_v(lUnit,Array,iStrt,iEnd,1,Ierr)
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
         iEnd  = iStrt + nPAM2Prim*nPAM2Bas - 1
         Do 288 iPrim = 0, nPAM2Prim - 1
            Call Read_v(lUnit,Array,iStrt+iPrim,iEnd,nPAM2Prim,ierr)
            If(ierr.ne.0) Then
               Call WarningMessage(2,
     &                         'GetBS: Error in reading GPA!!!')
               Call Abend()
            End if
 288     Continue
         iStrt = iEnd + 1
      End Do
*
      Call mma_allocate(dbsc(iCnttp)%PAM2,iEnd,Label='PAM2')
      dbsc(iCnttp)%PAM2(:)=Array(:)
*
      Call mma_deAllocate(Array)
      Return
      End
*
