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
* Copyright (C) 2003, Per-Olof Widmark                                 *
************************************************************************
************************************************************************
*                                                                      *
* This program tests the runfile utilities.                            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
* Written: July 2003                                                   *
*          Lund University, Sweden                                     *
*                                                                      *
************************************************************************
      Program TestRF
      Implicit None
#include "runtypes.fh"

      Integer      Mxdata
      Parameter    (MxData=64)

      Logical      UseOld

      Integer      iRc,iOpt
      Integer      iSeed
      Integer      Loop
      Integer      i

      Integer      nDataA,nDataB
      Integer      nDataC,nDataD
      Integer      nDataE,nDataF
      Integer      nDataAx,nDataBx
      Integer      nDataCx,nDataDx
      Integer      nDataEx,nDataFx
      Integer      RecTypA,RecTypB
      Integer      RecTypC,RecTypD
      Integer      RecTypE,RecTypF
      Real*8       DataA(MxData)
      Real*8       DataB(MxData)
      Integer      DataC(MxData)
      Integer      DataD(MxData)
      Character*1  DataE(MxData)
      Character*1  DataF(MxData)

      Real*8       Random_molcas
      External     Random_molcas

      Call NameRun('RUNFILE')
      UseOld=.false.
      If(UseOld) GoTo 500

      iSeed=12345
      nDataA=0
      nDataB=0
      nDataC=0
      nDataD=0
      nDataE=0
      nDataF=0

      Do i=1,MxData
         DataA(i)= 1.0d0*i
         DataB(i)=-1.0d0*i
         DataC(i)= i
         DataD(i)=-i
         DataE(i)= Char(64+mod(i,48))
         DataF(i)= Char(64+mod(i,48))
      End Do

      iRc=0
      iOpt=0
      Do Loop=1,16
         If(Random_molcas(iSeed).gt.0.75d0) Then
            nDataA=Int((MxData-4)*Random_molcas(iSeed)+4)
            Write(*,*) 'Adding A',nDataA
            Call dWrRun('Amat',DataA,nDataA)
         End If
         If(Random_molcas(iSeed).gt.0.75d0) Then
            nDataB=Int((MxData-4)*Random_molcas(iSeed)+4)
            Write(*,*) 'Adding B',nDataB
            Call dWrRun('Bmat',DataB,nDataB)
         End If
         If(Random_molcas(iSeed).gt.0.75d0) Then
            nDataC=Int((MxData-4)*Random_molcas(iSeed)+4)
            Write(*,*) 'Adding C',nDataC
            Call iWrRun('Cmat',DataC,nDataC)
         End If
         If(Random_molcas(iSeed).gt.0.75d0) Then
            nDataD=Int((MxData-4)*Random_molcas(iSeed)+4)
            Write(*,*) 'Adding D',nDataD
            Call iWrRun('Dmat',DataD,nDataD)
         End If
         If(Random_molcas(iSeed).gt.0.75d0) Then
            Call NameRun('RUNXXX')
            i=Int((MxData-4)*Random_molcas(iSeed)/4)
            nDataE=4*i+4
            Write(*,*) 'Adding E',nDataE
            Call cWrRun('Emat',DataE,nDataE)
            Call NameRun('RUNFILE')
         End If
         If(Random_molcas(iSeed).gt.0.75d0) Then
            i=Int((MxData-4)*Random_molcas(iSeed)/4)
            nDataF=4*i+4
            Write(*,*) 'Adding F',nDataF
            Call cWrRun('Fmat',DataF,nDataF)
         End If
      End Do

  500 Continue

      Do i=1,MxData
         DataA(i)= 0.0d0
         DataB(i)= 0.0d0
         DataC(i)= 0
         DataD(i)= 0
         DataE(i)= Char(0)
         DataF(i)= Char(0)
      End Do

      nDataAx=nDataA
      nDataBx=nDataB
      nDataCx=nDataC
      nDataDx=nDataD
      nDataEx=nDataE
      nDataFx=nDataF

      nDataA=0
      nDataB=0
      nDataC=0
      nDataD=0
      nDataE=0
      nDataF=0

      Call ffRun('aMat',nDataA,RecTypA)
      Call ffRun('bMat',nDataB,RecTypB)
      Call ffRun('cMat',nDataC,RecTypC)
      Call ffRun('dMat',nDataD,RecTypD)
      Call NameRun('RUNXXX')
      Call ffRun('eMat',nDataE,RecTypE)
      Call NameRun('RUNFILE')
      Call ffRun('fMat',nDataF,RecTypF)

      If(UseOld) Then
         nDataAx=nDataA
         nDataBx=nDataB
         nDataCx=nDataC
         nDataDx=nDataD
         nDataEx=nDataE
         nDataFx=nDataF
      End If

      If(nDataA.ne.nDataAx) Then
         Write(*,*) 'nDataA:',nDataA,nDataAx
         Stop
      End If

      If(nDataB.ne.nDataBx) Then
         Write(*,*) 'nDataB:',nDataB,nDataBx
         Stop
      End If

      If(nDataC.ne.nDataCx) Then
         Write(*,*) 'nDataC:',nDataC,nDataCx
         Stop
      End If

      If(nDataD.ne.nDataDx) Then
         Write(*,*) 'nDataD:',nDataD,nDataDx
         Stop
      End If

      If(nDataE.ne.nDataEx) Then
         Write(*,*) 'nDataE:',nDataE,nDataEx
         Stop
      End If

      If(nDataF.ne.nDataFx) Then
         Write(*,*) 'nDataF:',nDataF,nDataFx
         Stop
      End If

      If(RecTypA.ne.TypDbl) Then
         Write(*,*) 'RecTypA:',RecTypA,TypDbl
         Stop
      End If

      If(RecTypB.ne.TypDbl) Then
         Write(*,*) 'RecTypB:',RecTypB,TypDbl
         Stop
      End If

      If(RecTypC.ne.TypInt) Then
         Write(*,*) 'RecTypC:',RecTypC,TypInt
         Stop
      End If

      If(RecTypD.ne.TypInt) Then
         Write(*,*) 'RecTypD:',RecTypD,TypInt
         Stop
      End If

      If(RecTypE.ne.TypStr) Then
         Write(*,*) 'RecTypE:',RecTypE,TypStr
         Stop
      End If

      If(RecTypF.ne.TypStr) Then
         Write(*,*) 'RecTypF:',RecTypF,TypStr
         Stop
      End If

      If(nDataA.gt.0) Then
         Write(*,*) 'Testing A',nDataA
         Call dRdRun('amat',DataA,nDataA)
         Do i=1,nDataA
            If(Abs(DataA(i)-i).gt.1.0d-12) Then
               Write(*,*) 'A:',i,DataA(i)
            End If
         End Do
      End If

      If(nDataB.gt.0) Then
         Write(*,*) 'Testing B',nDataB
         Call dRdRun('bmat',DataB,nDataB)
         Do i=1,nDataB
            If(Abs(DataB(i)+i).gt.1.0d-12) Then
               Write(*,*) 'B:',i,DataB(i)
            End If
         End Do
      End If

      If(nDataC.gt.0) Then
         Write(*,*) 'Testing C',nDataC
         Call iRdRun('cmat',DataC,nDataC)
         Do i=1,nDataC
            If(DataC(i)-i.ne.0) Then
               Write(*,*) 'C:',i,DataC(i)
            End If
         End Do
      End If

      If(nDataD.gt.0) Then
         Write(*,*) 'Testing D',nDataD
         Call iRdRun('dmat',DataD,nDataD)
         Do i=1,nDataD
            If(DataD(i)+i.ne.0) Then
               Write(*,*) 'D:',i,DataD(i)
            End If
         End Do
      End If

      If(nDataE.gt.0) Then
         Call NameRun('RUNXXX')
         Write(*,*) 'Testing E',nDataE
         Call cRdRun('emat',DataE,nDataE)
         Do i=1,nDataE
            If(iChar(DataE(i))-64-mod(i,48).ne.0) Then
               Write(*,*) 'E:',i,DataE(i)
            End If
         End Do
         Call NameRun('RUNFILE')
      End If

      If(nDataF.gt.0) Then
         Write(*,*) 'Testing F',nDataF
         Call cRdRun('fmat',DataF,nDataF)
         Do i=1,nDataF
            If(iChar(DataF(i))-64-mod(i,48).ne.0) Then
               Write(*,*) 'F:',i,DataF(i)
            End If
         End Do
      End If

      Call DumpRun(iRc,iOpt)

      Stop
      End
