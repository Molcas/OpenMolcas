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
*                                                                      *
************************************************************************
      Subroutine SCF_Init()
************************************************************************
*                                                                      *
*     purpose: set up parameters that has to be predefined in SCF      *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
*
#include "mxdm.fh"
#include "file.fh"
#include "infscf.fh"
#include "infso.fh"
#include "llists.fh"
#include "twoswi.fh"
*
      Logical  Found, Reduce_Prt
      External Reduce_Prt
*
*----------------------------------------------------------------------*
*     Define files ( file names and unit numbers )                     *
*----------------------------------------------------------------------*
      LuOne = 10
      FnOne = 'ONEINT '
      LuOrd = 40
      FnOrd = 'ORDINT '
      LuOut = 20
      FnOut = 'SCFORB  '
      LuInp = 25
      LuDSt = 34
      FnDSt = 'DNSMAT  '
      LuOSt = 40
      FnOSt = 'DVXCDR  '
      LuTSt = 35
      FnTSt = 'TWOHAM  '
      LuGrd = 36
      FnGrd = 'GRADIENT'
      LuDGd = 37
      FnDGd = 'SODGRAD '
      Lux   = 38
      Fnx   = 'SOXVEC  '
      LuDel = 39
      FnDel = 'SODELTA '
      Luy   = 29
      Fny   = 'SOYVEC  '
*----------------------------------------------------------------------*
*     Initialize global variables                                      *
*----------------------------------------------------------------------*
      Call Peek_iScalar('nSym',nSym)
      nAtoms = 0
      Call ICopy(MxSym,[0],0,nBas,1)
      Call ICopy(MxSym,[0],0,nOrb,1)
      Call ICopy(MxSym*2,[0],0,nOcc,1)
      If(iUHF.eq.0) Then
         Call Put_iArray('nIsh',nOcc(1,1),nSym)
      Else
         Call Put_iArray('nIsh',nOcc(1,1),nSym)
         Call Put_iArray('nIsh beta',nOcc(1,2),nSym)
      Endif
      Call ICopy(MxSym,[0],0,nFro,1)
      Call Put_iArray('nFro',nFro,nSym)
      Call ICopy(MxSym,[0],0,nFrz,1)
      Call ICopy(MxSym,[0],0,nDel,1)
      Call qpg_iarray('nDel',Found,ndata)
      If (.not.Found) Then
         Call Put_iArray('nDel',nDel,nSym)
      End If
      Thize = 1.0d-6
      EThr  = 1.0d-9
C     DThr  = 1.0d-5
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
      DltNTh= 0.2d-4
C     FThr   =  0.5d-6
      FThr   =  1.5d-4
      QudThr = 1.0d-5
      Energy = 0.0D0
      E1     = 0.0D0
      E2     = 0.0D0
      EKin   = 0.0D0
      PotNuc = 0.0D0
      EneV = 0.0d0
      E1V  = 0.0d0
      E2V  = 0.0d0
      EKinV= 0.0d0
      OnlyProp=.false.
      NoProp=.false.
      FckAuf=.True.
*
* Order: SCF(0), Guessorb(1), Lumorb(2), Core(4)
* New order: SCF(0), Lumorb(2), Guessorb(1), Core(4)
*
      InVec   = 1
      Do i=1,nStOpt
         LstVec(i)=-1
      End Do
      LstVec(1)=0
      LstVec(2)=2
      LstVec(3)=1
      LstVec(4)=4
      nIter(0)= MxIter
      nIter(1)= MxIter
      nIterP  = 1
      iter    = 1
      iter0   = 0
      iterso  = 0
*
      iPrint=iPrintLevel(-1)
      jPrint=iPrint
      If (Reduce_Prt().and.iPrint.lt.3) jPrint=0
*
      iPrOrb = 1
      kIvo   = 0
      iUHF   = 0
      iROHF  = 0
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
      Call ICopy(MxOptm,[-1],0,kDisk,1)
      Call ICopy(MxDDsk,[-1],0,iDisk,1)
      Call ICopy(MxKeep,[0],0,MapDns,1)
      FrstDs = .True.
      FrstDa = .True.
      PreSch = .False.
      MiniDn = .True.
      WrOutD = .False.
      c1Diis = .False.
      RSRFO  = .False.
      Scrmbl = .False.
      RFpert = .False.
      PmTime = .False.
      EmConv = .False.
      AddFragments = .False.
      call dcopy_(nFld,[0.0d0],0,TimFld,1)
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
      ExFac=1.0D0
      WarnCfg=.False.
      WarnPocc=.False.
      WarnSlow=.False.
*
      Call Get_iScalar('FMM',iFMM)
      DoFMM=iFMM.eq.1
*
*     Initialize energy due to external potential on the run file. This
*     to make sure that it can be read unconditionally and is zero if
*     not redefined by DrvXV
      Call Poke_dScalar('KSDFT energy',0.0D0)
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
      End
