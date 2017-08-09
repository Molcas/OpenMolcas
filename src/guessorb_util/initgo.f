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
* Copyright (C) 2004, Per-Olof Widmark                                 *
************************************************************************
************************************************************************
*                                                                      *
* This routine initialize guessorb.                                    *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University                                             *
*          Sweden                                                      *
* Written: Oct 2004                                                    *
*                                                                      *
************************************************************************
      Subroutine InitGO(StandAlone)
      Implicit None
      Logical StandAlone
#include "Molcas.fh"
#include "commgo.fh"
*----------------------------------------------------------------------*
* Local variables                                                      *
*----------------------------------------------------------------------*
      Logical Debug
      Logical Trace
      Logical Found
      Integer nBasTot
      Integer iSym
      Integer iBas
      Integer iPrt
      Integer i
*----------------------------------------------------------------------*
* External entry points                                                *
*----------------------------------------------------------------------*
      Integer  iPrintLevel
      External iPrintLevel
*----------------------------------------------------------------------*
* Setup                                                                *
*----------------------------------------------------------------------*
      Debug=.false.
      Trace=.false.
      If(Trace) Write(6,*) '>>> Entering initgo'
*----------------------------------------------------------------------*
* Set default for MO printing.                                         *
*----------------------------------------------------------------------*
      iPrt=iPrintLevel(-1)
      If(iPrt.ge.4) Then
         PrintMOs=.true.
         PrintEor=.true.
         PrintPop=.true.
         iPrFmt=3
         PrThr=1.0d6
      Else If(iPrt.ge.3) Then
         PrintMOs=.false.
         PrintEor=.false.
         PrintPop=.false.
         iPrFmt=1
         PrThr=5.0d0
      Else
         PrintMOs=.false.
         PrintEor=.false.
         PrintPop=.false.
         PrThr=5.0d0
      End If
*----------------------------------------------------------------------*
* Set other defaults.                                                  *
*----------------------------------------------------------------------*
      write(6,*) "initgo1 yes"
      Call Qpg_dScalar('S delete thr',Found)
      write(6,*) "initgo2 yes"
      If(Found) Then
      write(6,*) "initgo3"
         Call Get_dScalar('S delete thr',SThr)
      write(6,*) "initgo4"
      Else
      write(6,*) "initgo5 yes"
         SThr=1.0d-9
         Call Put_dScalar('S delete thr',SThr)
      write(6,*) "initgo6 yes"
      End If
      write(6,*) "initgo7 yes"
      Call Qpg_dScalar('T delete thr',Found)
      write(6,*) "initgo8 yes"
      If(Found) Then
      write(6,*) "initgo9"
         Call Get_dScalar('T delete thr',TThr)
      write(6,*) "initgo10"
      Else
         TThr=1.0d+6
      write(6,*) "initgo11 yes"
         Call Put_dScalar('T delete thr',TThr)
      write(6,*) "initgo12 yes"
      End If
      GapThr=0.01d0
*----------------------------------------------------------------------*
* Get basic data from runfile.                                         *
*----------------------------------------------------------------------*
      write(6,*) "initgo13 yes"
      Call get_iscalar('nSym',nSym)
      write(6,*) "initgo14 yes"
      Call get_iarray('nBas',nBas,nSym)
      write(6,*) "initgo15 yes"
      Do iSym=1,MxSym
         nOcc(iSym)=0
         nVir(iSym)=0
         nDel(iSym)=0
      End Do
      nBasTot=0
      Do iSym=1,nSym
         nBasTot=nBasTot+nBas(iSym)
      End Do
*VB      If(Debug) Then
         Write(6,'(a,8i5)') 'initgo: nSym',nSym
         Write(6,'(a,8i5)') 'initgo: nBas',nBas
         Write(6,'(a,8i5)') 'initgo: nOcc',nOcc
         Write(6,'(a,8i5)') 'initgo: nVir',nVir
         Write(6,'(a,8i5)') 'initgo: nBasTot',nBasTot
*VB      End If
      write(6,*) "initgo16 yes"
      Call get_iscalar('Unique Atoms',nNuc)
      write(6,*) "initgo17 yes"
      If(nNuc.gt.MxAtom) Then
      write(6,*) "initgo18"
         Call SysAbendMsg('initgo','Fatal:',
     &                    'Too many atoms, increase MxAtom')
      End If
      write(6,*) "initgo19 yes"
      Call get_carray('Unique Atom Names',Name,LENIN*nNuc)
      write(6,*) "initgo20 yes"
      Call get_carray('Unique Basis Names',Label,(LENIN4)*nBasTot)
      write(6,*) "initgo21 yes"
      Call get_darray('Nuclear Charge',xCharge,nNuc)
      write(6,*) "initgo22 yes"
*VB      If(Debug) Then
         Write(6,'(a,8i5)')    'initgo: nNuc',nNuc
         Write(6,'(a,8i5)')    'initgo: nBasTot',nBasTot
         Write(6,'(a,8f12.6)') 'initgo: Charge',(xCharge(i),i=1,nNuc)
         Write(6,'(a,8a4)')    'initgo: Name ',(Name(i),i=1,nNuc)
         Write(6,'(a)') 'initgo: Basis functions'
         Do iBas=1,nBasTot
            Write(6,'(2a)') Label(iBas)(1:LENIN),
     &                      Label(iBas)(LENIN1:LENIN4)
         End Do
*VB      End If
*----------------------------------------------------------------------*
* Done                                                                 *
*----------------------------------------------------------------------*
      If(Trace) Write(6,*) '<<< Exiting initgo'
      Return
c Avoid unused argument warnings
      write(6,*) "initgo23"
      If (.False.) Call Unused_logical(StandAlone)
      write(6,*) "initgo24"
      End
