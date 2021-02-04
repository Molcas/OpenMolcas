!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2004, Per-Olof Widmark                                 *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This routine initialize guessorb.                                    *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: Oct 2004                                                    *
!                                                                      *
!***********************************************************************
      Subroutine InitGO(StandAlone)
      use GuessOrb_Global, only: GapThr, iPrFmt, Label, LenIn, LenIn1, LenIn8, MxAtom, MxSym, Name, nBas, nDel, nNuc, nOcc, nSym, &
                                 nVir, PrintEor, PrintMOs, PrintPop, PrThr, SThr, TThr, xCharge
      Implicit None
      Logical StandAlone
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
      Logical Debug
      Logical Trace
      Logical Found
      Integer nBasTot
      Integer iSym
      Integer iBas
      Integer iPrt
      Integer i
!----------------------------------------------------------------------*
! External entry points                                                *
!----------------------------------------------------------------------*
      Integer  iPrintLevel
      External iPrintLevel
!----------------------------------------------------------------------*
! Setup                                                                *
!----------------------------------------------------------------------*
      Debug=.false.
      Trace=.false.
      If(Trace) Write(6,*) '>>> Entering initgo'
!----------------------------------------------------------------------*
! Set default for MO printing.                                         *
!----------------------------------------------------------------------*
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
!----------------------------------------------------------------------*
! Set other defaults.                                                  *
!----------------------------------------------------------------------*
      Call Qpg_dScalar('S delete thr',Found)
      If(Found) Then
         Call Get_dScalar('S delete thr',SThr)
      Else
         SThr=1.0d-9
         Call Put_dScalar('S delete thr',SThr)
      End If
      Call Qpg_dScalar('T delete thr',Found)
      If(Found) Then
         Call Get_dScalar('T delete thr',TThr)
      Else
         TThr=1.0d+6
         Call Put_dScalar('T delete thr',TThr)
      End If
      GapThr=0.01d0
!----------------------------------------------------------------------*
! Get basic data from runfile.                                         *
!----------------------------------------------------------------------*
      Call get_iscalar('nSym',nSym)
      Call get_iarray('nBas',nBas,nSym)
      Do iSym=1,MxSym
         nOcc(iSym)=0
         nVir(iSym)=0
         nDel(iSym)=0
      End Do
      nBasTot=0
      Do iSym=1,nSym
         nBasTot=nBasTot+nBas(iSym)
      End Do
      If(Debug) Then
         Write(6,'(a,8i5)') 'initgo: nSym',nSym
         Write(6,'(a,8i5)') 'initgo: nBas',nBas
         Write(6,'(a,8i5)') 'initgo: nOcc',nOcc
         Write(6,'(a,8i5)') 'initgo: nVir',nVir
         Write(6,'(a,8i5)') 'initgo: nBasTot',nBasTot
      End If
      Call get_iscalar('Unique Atoms',nNuc)
      If(nNuc.gt.MxAtom) Then
         Call SysAbendMsg('initgo','Fatal:',                            &
     &                    'Too many atoms, increase MxAtom')
      End If
      Call get_carray('Unique Atom Names',Name,LENIN*nNuc)
      Call get_carray('Unique Basis Names',Label,(LENIN8)*nBasTot)
      Call get_darray('Nuclear Charge',xCharge,nNuc)
      If(Debug) Then
         Write(6,'(a,8i5)')    'initgo: nNuc',nNuc
         Write(6,'(a,8i5)')    'initgo: nBasTot',nBasTot
         Write(6,'(a,8f12.6)') 'initgo: Charge',(xCharge(i),i=1,nNuc)
         Write(6,'(a,8a4)')    'initgo: Name ',(Name(i),i=1,nNuc)
         Write(6,'(a)') 'initgo: Basis functions'
         Do iBas=1,nBasTot
            Write(6,'(2a)') Label(iBas)(1:LENIN),                       &
     &                      Label(iBas)(LENIN1:LENIN8)
         End Do
      End If
!----------------------------------------------------------------------*
! Done                                                                 *
!----------------------------------------------------------------------*
      If(Trace) Write(6,*) '<<< Exiting initgo'
      Return
! Avoid unused argument warnings
      If (.False.) Call Unused_logical(StandAlone)
      End
