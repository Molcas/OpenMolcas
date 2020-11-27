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
      SubRoutine Status(kIter,Energy,rGrad,GrdMax,GrdLbl,StpMax,
     &                  StpLbl,Ex,nLines,delE,iNeg,UpMeth,HUpMet,
     &                  Step_Trunc,Print_Status)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
*     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
*             University of Lund, SWEDEN                               *
*             May '91                                                  *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
#include "stdalloc.fh"
      Character*128 GrdLbl*8, StpLbl*8, UpMeth*6, HUpMet*8, Step_Trunc*1
      Character*8 lNeg
      Integer iNeg(2)
      Logical Print_Status
      Character(LEN=128), Allocatable:: Lines(:)
*                                                                      *
************************************************************************
*                                                                      *
      Interface
      Subroutine Put_cArray(Label,Data,nData)
      Character*(*) Label
      Integer       nData
      Character*(*) Data
      End Subroutine Put_cArray
      End Interface
*                                                                      *
************************************************************************
*                                                                      *

      Lu=6
      iRout = 52
      iPrint = nPrint(iRout)
*
      Call mma_allocate(Lines,[-1,nLines],Label='Lines')
*
*     Pick up previous energy
*
      If (kIter.eq.1) Then
         Lines(:)=' '
         iter=1
         Write (Lines(-1),'(A)')
     &    '                       Energy '//
     &    '    Grad      Grad          '//
     &    '    Step                 Estimated   Geom'//
     &    '       Hessian'
         Write (Lines(0),'(A)')
     &    'Iter      Energy       Change '//
     &    '    Norm      Max    Element '//
     &    '   Max     Element     Final Energy Update'//
     &    ' Update   Index'
      Else
         Call Get_cArray('Slapaf Info 3',Lines,(nLines+2)*128)
*--------Find first blank line
         iter = kIter
      End If
*
      If (iter.gt.nLines) Then
         Call WarningMessage(2,'Status: iter.gt.nLines')
         Write (Lu,*) 'iter,nLines=',iter,nLines
         Call abend()
      Else If (iter.lt.1) Then
         Call WarningMessage(2,'Status: iter.lt.1')
         Call abend()
      End If
*
      Write(lNeg,'(I3)') iNeg(1)
      If (iNeg(2).ne.iNeg(1)) Then
         If (iNeg(2).gt.99) Then
            Write(lNeg(4:8),'("(",I3,")")') iNeg(2)
         Else If (iNeg(2).gt.9) Then
            Write(lNeg(4:7),'("(",I2,")")') iNeg(2)
         Else
            Write(lNeg(4:6),'("(",I1,")")') iNeg(2)
         End If
      End If
      Write (Lines(iter),
     &    '(I3,F16.8,F12.8,2(F9.6,1X),A8,F9.6,A1,1X,A8,F16.8,1X,'//
     &    'A,1X,A,1X,A)')
     &    iter,Energy,delE,rGrad,GrdMax,GrdLbl,StpMax,Step_Trunc,
     &    StpLbl,Ex,UpMeth,HUpMet,lNeg
*                                                                      *
************************************************************************
*                                                                      *
      Lu_file=isfreeunit(8)
*
*     Turn off updating the structure file during numerical
*     integrations steps.
*
      nvv = 1
      If (Print_Status) nvv=2
      Do ivv=1,nvv
         If (ivv.eq.1) then
           Lu_out=Lu
         Else
            Call Molcas_Open(Lu_File,'STRUCTURE')
           Lu_out=Lu_file
         End If
*
         If (ivv.eq.2.or.iPrint.ge.5) Then
            Write (Lu_out,*)
            Write (Lu_out,'(A)') '***********************'//
     &         '*******************************************************'
     &         //'****************************************'
            Write (Lu_out,'(A)') '*                              '//
     &         '      Energy Statistics for Geometry Optimization      '
     &                    //'                               *'
            Write (Lu_out,'(A)') '***********************'//
     &         '*******************************************************'
     &         //'****************************************'
            Do i = -1, iter
               Write (Lu_out,'(A118)') Lines(i)
            End Do
            Write (Lu_out,*)
         End If
         If (ivv.eq.2) Close(Lu_file)
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Call Put_cArray('Slapaf Info 3',Lines(-1),(nLines+2)*128)
      Call mma_deallocate(Lines)
      Return
      End
