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
!                                                                      *
!***********************************************************************
      Subroutine SCF_Init()
!***********************************************************************
!                                                                      *
!     purpose: set up parameters that has to be predefined in SCF      *
!                                                                      *
!***********************************************************************
      Use InfSO, only: qNRTh, DltnTH, Energy, IterSO
      use InfSCF, only: nAtoms, nOcc, nD, Thize, EThr, DThr, DelThr, DIISTh, FThr, QudThr, E1, E2, EKin, OnlyProp,     &
                        NoProp, FckAuf, InVec, nIterP, Iter, iPrint, jPrint, iPrOrb, kIVO, iCoCo, jVOut,                 &
                        DIIS, Damping, One_Grid, Two_Thresholds, IDKeep, lPaper, nDens, kOptim, AccCon, nDisc,           &
                        nCore, kDisk, PotNuc, EneV, E1V, E2V, iDisk, EkInv, MapDns, PreSch, MiniDn, WrOutD, C1DIIS,      &
                        RSRFO, RGEK, Scrmbl, RFPert, pmTime, EmConv, AddFragments, rTemp, TemFac, TStop, KSDFT,          &
                        ExFac, WarnCFG, WarnPOCC, WarnSlow, DoFMM, nIter, NamFld, LstVec, nBas, nDel, nFro,              &
                        nOrb, nSym, TimFld, nFrz
      use Constants, only: Zero, One
      use MxDM, only: MxIter
      use Constants, only: Zero, One
      Implicit None
!
#include "twoswi.fh"
#include "hfc_logical.fh"
!
      Integer  iFMM, iPrintLevel, nData
      Logical  Found, Reduce_Prt
      External Reduce_Prt
!----------------------------------------------------------------------*
!     Initialize global variables                                      *
!----------------------------------------------------------------------*
      Call Peek_iScalar('nSym',nSym)
      nAtoms = 0
      nBas(:)=0
      nOrb(:)=0
      nOcc(:,:)=0
      If(nD==1) Then
         Call Put_iArray('nIsh',nOcc(1,1),nSym)
      Else
         Call Put_iArray('nIsh',nOcc(1,1),nSym)
         Call Put_iArray('nIsh beta',nOcc(1,2),nSym)
      Endif
      nFro(:)=0
      Call Put_iArray('nFro',nFro,nSym)
      nFrz(:)=0
      nDel(:)=0
      Call qpg_iarray('nDel',Found,ndata)
      If (.not.Found) Then
         Call Put_iArray('nDel',nDel,nSym)
      End If
      Thize = 1.0d-6
      EThr  = 1.0d-9
      DThr  = 1.0d-4
      Call Qpg_dScalar('S delete thr',Found)
      If(Found) Then
         Call Get_dScalar('S delete thr',DelThr)
      Else
         DelThr= 1.0d-5
         Call Put_dScalar('S delete thr',DelThr)
      End If
      DiisTh= 0.15d+00
      QNRTh = 0.075d+00
!     DltNTh= 0.2d-4
      DltNTh= 0.1d-2
      FThr   =  1.5d-4
      QudThr = 1.0d-5
      Energy = Zero
      E1     = Zero
      E2     = Zero
      EKin   = Zero
      PotNuc = Zero
      EneV = Zero
      E1V  = Zero
      E2V  = Zero
      EKinV= Zero
      OnlyProp=.false.
      NoProp=.false.
      FckAuf=.True.
!
! Order: SCF(0), Guessorb(1), Lumorb(2), Core(4)
! New order: SCF(0), Lumorb(2), Guessorb(1), Core(4)
!
      InVec   = 1
      LstVec(:)=-1
      LstVec(1)=0
      LstVec(2)=2
      LstVec(3)=1
      LstVec(4)=4
      nIter(0)= MxIter
      nIter(1)= MxIter
      nIterP  = 1
      iter    = 1
      iterso  = 0
!
      iPrint=iPrintLevel(-1)
      jPrint=iPrint
      If (Reduce_Prt().and.iPrint.lt.3) jPrint=0
!
      iPrOrb = 1
      kIvo   = 0
      nD = 1
      UHF_HFC   = .False.
      iCoCo = 0
      jVOut = 2
      Diis = .True.
      Damping = .True.
      One_Grid=.False.
      Two_Thresholds=.True.
      iDKeep = 4
      lPaper = 132
      nDens = 0
      kOptim = 1
      AccCon = '         '
      nDisc  =  2000
      nCore = 512
      kDisk(:)=-1
      iDisk(:,1)=-1
      MapDns(:)=0
      PreSch = .False.
      MiniDn = .True.
      WrOutD = .False.
      c1Diis = .False.
      RSRFO  = .False.
      RGEK   = .False.
      Scrmbl = .False.
      RFpert = .False.
      PmTime = .False.
      EmConv = .False.
      AddFragments = .False.
      TimFld(:)=Zero
      NamFld( 1)='1) Input processing                         :'
      NamFld( 2)='2) Wave function optimization               :'
      NamFld( 3)='     Line Search Iterations    (QNR steps)  :'
      NamFld( 4)='  a ) calculation of the density            :'
      NamFld( 5)='  b ) contraction with integrals            :'
      NamFld( 6)='  c ) acceleration of convergence           :'
      NamFld( 7)='        recursive BFGS         (QNR steps)  :'
      NamFld( 8)='  d ) solution to Roothaan-Hall equations   :'
      NamFld( 9)='  d'') rotate MOs C with U      (QNR steps)  :'
      NamFld(10)='        U=exp(kap)                          :'
      NamFld(11)='          via Taylor expansion (sin/cos)    :'
      NamFld(12)='          via transformation to Schur basis :'
      NamFld(13)='  e'') transf. Fck Mat. with C  (QNR steps)  :'
      NamFld(14)='  f ) other calculations                    :'
      NamFld(15)='3) Final processing (generation of outputs) :'
      NamFld(16)='   T O T A L                                :'
      NDDO = .False.
      RTemp = 0.04D0
      TemFac= 0.02D0
      TStop = 0.00D0
      KSDFT='SCF '
      ExFac=One
      WarnCfg=.False.
      WarnPocc=.False.
      WarnSlow=.False.
!
      Call Get_iScalar('FMM',iFMM)
      DoFMM=iFMM.eq.1
!
!     Initialize energy due to external potential on the run file. This
!     to make sure that it can be read unconditionally and is zero if
!     not redefined by DrvXV
      Call Poke_dScalar('KSDFT energy',Zero)
!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
      End Subroutine SCF_Init
