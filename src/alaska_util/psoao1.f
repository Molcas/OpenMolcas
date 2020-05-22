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
* Copyright (C) 1990,1992, Roland Lindh                                *
*               1990, IBM                                              *
************************************************************************
      SubRoutine PSOAO1(nSO,MemPrm,MemMax,
     &                            iAnga, iCmpa, iShela, iFnc,
     &                            iBas,  iBsInc, jBas,  jBsInc,
     &                            kBas,  kBsInc, lBas,  lBsInc,
     &                            iPrim, iPrInc, jPrim, jPrInc,
     &                            kPrim, kPrInc, lPrim, lPrInc,
     &                            ipMem1,ipMem2,
     &                            Mem1,  Mem2,  MemPSO)
************************************************************************
*                                                                      *
*  Object: to partion the SO and AO block. It will go to some length   *
*          before it will start and break up the SO block. This will   *
*          reduce the total flop count. However, as we are breaking up *
*          the AO block this will affect the vectorization. Hence, at  *
*          some point it will actually be better to recompute the      *
*          primitives.                                                 *
*          Current stratergy:                                          *
*          1. Start reducing the length of the primitives in the order *
*             lPrim,jPrim.                                             *
*          2. Reduce the size of the SO block by reducing the number of*
*             basis functions in the order lBas, jBas.                 *
*          3. Terminate run telling job max and min of additional      *
*             memory needed to perform the calculation.                *
*                                                                      *
* Called from: Drvg1                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              Change                                                  *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*                                                                      *
*             Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             Modified to first order derivatives. January '92         *
************************************************************************
      use aces_stuff, only: nGamma, Gamma_On
      use PSO_Stuff
      Implicit Real*8 (A-H,O-Z)
#include "WrkSpc.fh"
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "lCache.fh"
#include "pstat.fh"
      Integer   iAnga(4), iCmpa(4), nPam(4,0:7), iiBas(4),
     &          iShela(4), iFnc(4)
      Logical QiBas, QjBas, QkBas, QlBas, QjPrim, QlPrim, Fail
      Integer   iTwoj(0:7)
      Data iTwoj/1,2,4,8,16,32,64,128/
*
*     Statement function to compute canonical index
*
      nElem(i) = (i+1)*(i+2)/2
*
      iRout = 10
      iPrint = nPrint(iRout)
*     Call qEnter('PSOAO1')
      la = iAnga(1)
      lb = iAnga(2)
      lc = iAnga(3)
      ld = iAnga(4)
      iCmp = iCmpa(1)
      jCmp = iCmpa(2)
      kCmp = iCmpa(3)
      lCmp = iCmpa(4)
      iTotal = iTotal + 1
      mabcd=nElem(la)*nElem(lb)*nElem(lc)*nElem(ld)
      nabcd=iCmp*jCmp*kCmp*lCmp
*
      If (force_part_c) Then
         iBsInc = (iBas+1)/2
         jBsInc = (jBas+1)/2
         kBsInc = (kBas+1)/2
         lBsInc = (lBas+1)/2
      Else
         iBsInc = iBas
         jBsInc = jBas
         kBsInc = kBas
         lBsInc = lBas
      End If
      If (force_part_p) Then
         jPrInc = (jPrim+1)/2
         lPrInc = (lPrim+1)/2
      Else
         jPrInc = jPrim
         lPrInc = lPrim
      End If
      iPrInc = iPrim
      kPrInc = kPrim
*
*     Call GetMem(' ','MAX ','REAL',iDum,MemMax)
 999  Continue
      QjPrim = .False.
      QlPrim = .True.
      QiBas  = .False.
      QjBas  = .False.
      QkBas  = .False.
      QlBas  = .False.
      Mem0 = MemMax
*
*-----*** Work1 ***
*
*-----Memory for 2nd order density matrix in SO basis.
*
      kSOInt = nSO*iBsInc*jBsInc*kBsInc*lBsInc
      Mem1 = kSOInt
*
*-----Allocate memory for MO to SO/AO transformation
*     of the 2nd order density matrix for this shell quadruplet.
*
      If (lPSO) Then
         iiBas(1) = iBsInc
         iiBas(2) = jBsInc
         iiBas(3) = kBsInc
         iiBas(4) = lBsInc
         Call ICopy(4*8,[0],0,nPam,1)
         MemPSO = 1
         nTmp2 = 0
*
         Do jPam = 1, 4
            iTmp1= 0
            nTmp1= 0
            Do j = 0, nIrrep-1
               Do i1 = 1, iCmpa(jPam)
                  If (iAnd(IrrCmp(IndS(iShela(jPam))+i1),
     &                iTwoj(j)).ne.0) Then
                      nPam(jPam,j) = nPam(jPam,j) + iiBas(jPam)
                      nTmp1= nTmp1+ iiBas(jPam)
                      iTmp1= iTmp1+ 1
                  End If
               End Do
            End Do
            MemPSO = MemPSO * nTmp1
            nTmp2 = nTmp2 + nTmp1
            iFnc(jPam) = iTmp1
         End Do
         MemScr=MemTra(nPam)
         nFac = 4
         nTmp2 = nTmp2 + 4
      Else
         MemScr=0
         MemPSO=0
         nFac = 0
         nTmp2 = 0
      End If
      MemAux0= MemPSO + MemScr + nFac*nDim + nTmp2 + 4
      If (Mem1+1+MemAux0.gt.Mem0) Then
         MaxReq=Max(MaxReq,Mem1+1+MemAux0-Mem0)
         QjPrim = .False.
         QlPrim = .False.
         QiBas  = .False.
         QjBas  = .False.
         QkBas  = .False.
         QlBas  = .True.
         Call Change(iBas, iBsInc,QiBas, kBas, kBsInc,QkBas,
     &               jBas, jBsInc,QjBas, lBas, lBsInc,QlBas,
     &               jPrim,jPrInc,QjPrim,lPrim,lPrInc,QlPrim,MaxReq,
     &               Fail)
         If (Fail) Then
            Write (6,*) ' Memory allocation failed for Work1'
            Write (6,'(2I3,L1,2I3,L1)')
     &            iBas, iBsInc,QiBas, kBas, kBsInc,QkBas
            Write (6,'(2I3,L1,2I3,L1)')
     &            jBas, jBsInc,QjBas, lBas, lBsInc,QlBas
            Write (6,'(2I3,L1,2I3,L1)')
     &            jPrim,jPrInc,QjPrim,lPrim,lPrInc,QlPrim
            Write (6,*) MemMax,Mem0,Mem1,MemAux0+1
            Write (6,*) MemPSO,MemScr,4*nDim,nTmp2+4
            Call Abend()
         End If
         Go To 999
      End If
*-----Subtracte one additional word for getmem's internal error check
      Mem0 = Mem0 - Mem1 - 1
*
*-----*** Work2 and Work4 ***
*
*-----Memory for 2nd order density matrix in contracted basis (both
*     cartesian and spherical harmonic) and in primitive basis.
*     MemDeP: Target for desymmetrization
*     MemTrn: Scratch and target for decontraction
*     MemAux: Contracted 2nd order density matrix (if partial decon.)
*     MemSph: transformation spherical harmonics to cartesian, source
*             and target.
*     nGamma: temporary storage of bins as read from file.
*
      If (Gamma_On) Then
         iiBas(1) = iBas
         iiBas(2) = jBas
         iiBas(3) = kBas
         iiBas(4) = lBas
         nGamma = 1
*
         Do jPam = 1, 4
            nTmp1= 0
            Do j = 0, nIrrep-1
               Do i1 = 1, iCmpa(jPam)
                  If (iAnd(IrrCmp(IndS(iShela(jPam))+i1),
     &                iTwoj(j)).ne.0) Then
                      nTmp1= nTmp1+ iiBas(jPam)
                  End If
               End Do
            End Do
            nGamma = nGamma * nTmp1
         End Do
      Else
         nGamma=0
      End If
*
      MemDeP = nabcd * iBsInc*jBsInc*kBsInc*lBsInc
*
      MemTrn = mabcd * Max(iBsInc*jBsInc*kBsInc*lBsInc,
     &                     iPrInc*jPrInc*kPrInc*lPrInc)
*-----If partial decontraction we need to keep the contracted 2nd
*     order density matrix. (Work4)
      If (jPrInc.ne.jPrim .or. lPrInc.ne.lPrim) Then
         MemAux = mabcd*iBsInc*jBsInc*kBsInc*lBsInc
      Else
         MemAux = 0
      End If
      MemSph = mabcd * iBsInc*jBsInc*kBsInc*lBsInc
      Mem2 = Max(MemTrn+MemAux,MemDeP,MemSph,nGamma+MemAux0)
      If (Mem2+1.gt.Mem0) Then
         MaxReq=Max(MaxReq,Mem2+1-Mem0)
         Call Change(iBas, iBsInc,QiBas, kBas, kBsInc,QkBas,
     &               jBas, jBsInc,QjBas, lBas, lBsInc,QlBas,
     &               jPrim,jPrInc,QjPrim,lPrim,lPrInc,QlPrim,MaxReq,
     &               Fail)
         If (Fail) Then
            Write (6,*) ' Memory allocation failed for Work2'
            Write (6,'(2I3,L1,2I3,L1)')
     &            iBas, iBsInc,QiBas, kBas, kBsInc,QkBas
            Write (6,'(2I3,L1,2I3,L1)')
     &            jBas, jBsInc,QjBas, lBas, lBsInc,QlBas
            Write (6,'(2I3,L1,2I3,L1)')
     &            jPrim,jPrInc,QjPrim,lPrim,lPrInc,QlPrim
            Write (6,*) MemMax,Mem0,Mem1,Mem2
            Call Abend()
         End If
         Go To 999
      End If
*-----Subtracte one additional word for getmem's internal error check
      Mem0 = Mem0 - Mem2 - 1
*
*-----*** Work3 and Work5 ***
*
*-----Scratch for decontraction and transformation to spherical gaussian.
*     Working array for Rysg1.
*     Scratch area for resolving degeneracies due to the double coset
*     treatement of the symmetry.
*     MemTrn: Scratch for decontraction
*     MemRys: Scratch for calualation of primitive integral gradients.
*     MemScr: Scratch for prescreening
*
*     If mabcd=/=1 we need to do a transpose!
      iFac = 0
      If (mabcd.ne.1) iFac = 1
*
      nCache = (3*lCache)/4 - iBas*iPrim - jBas*jPrim
      lSize= iBsInc*jBsInc + Max(iBsInc*jPrInc,iPrInc*jBsInc)
      nVec1 = kBsInc*lBsInc * mabcd
      IncVec = Min(Max(1,nCache/lSize),nVec1)
      nA2 = Max(iBsInc*jPrInc,iPrInc*jBsInc)*IncVec
      nA3 = iPrInc*jPrInc*nVec1
      MemTrn=nA2+nA3
*
      nCache = (3*lCache)/4 - kBas*kPrim - lBas*lPrim
      lSize= kBsInc*lBsInc + Max(kBsInc*lPrInc,kPrInc*lBsInc)
      nVec1 = iPrInc*jPrInc * mabcd
      IncVec = Min(Max(1,nCache/lSize),nVec1)
      nA2 = Max(kBsInc*lPrInc,kPrInc*lBsInc)*IncVec
      MemTrn=Max(MemTrn,nA3+nA2,
     &           iFac*mabcd*iBsInc*jBsInc*kBsInc*lBsInc)
*
************************************************************************
*
      nCache = (3*lCache)/4 - kBas*kPrim - lBas*lPrim
      lSize= kBsInc*lBsInc + Max(kBsInc*lPrInc,kPrInc*lBsInc)
      nVec1 = iBsInc*jBsInc * mabcd
      IncVec = Min(Max(1,nCache/lSize),nVec1)
      nA2 = Max(kBsInc*lPrInc,kPrInc*lBsInc)*IncVec
      nA3 = kPrInc*lPrInc*nVec1
      MemTrn=Max(MemTrn,nA3+nA2)
*
      MemTrn=Max(MemTrn,nA2+nA3)
      nCache = (3*lCache)/4 - iBas*iPrim - jBas*jPrim
      lSize= iBsInc*jBsInc + Max(iBsInc*jPrInc,iPrInc*jBsInc)
      nVec1 = kPrInc*lPrInc * mabcd
      IncVec = Min(Max(1,nCache/lSize),nVec1)
      nA2 = Max(iBsInc*jPrInc,iPrInc*jBsInc)*IncVec
*
      MemTrn=Max(MemTrn,nA2+nA3,
     &           iFac*mabcd*iBsInc*jBsInc*kBsInc*lBsInc)
*
      MemRys=MemPrm * iPrInc*jPrInc*kPrInc*lPrInc
      MemScr=(2*mabcd+1)*iPrInc*jPrInc*kPrInc*lPrInc
     &      +iPrInc*jPrInc+kPrInc*lPrInc
      Mem3 = Max(MemTrn, MemRys, MemScr)
      If (Mem3+1.gt.Mem0) Then
         MaxReq=Max(MaxReq,Mem3+1-Mem0)
         Call Change(iBas, iBsInc,QiBas, kBas, kBsInc,QkBas,
     &               jBas, jBsInc,QjBas, lBas, lBsInc,QlBas,
     &               jPrim,jPrInc,QjPrim,lPrim,lPrInc,QlPrim,MaxReq,
     &               Fail)
         If (Fail) Then
            Write (6,*) ' Memory allocation failed for Work3'
            Write (6,'(2I3,L1,2I3,L1)')
     &            iBas, iBsInc,QiBas, kBas, kBsInc,QkBas
            Write (6,'(2I3,L1,2I3,L1)')
     &            jBas, jBsInc,QjBas, lBas, lBsInc,QlBas
            Write (6,'(2I3,L1,2I3,L1)')
     &            jPrim,jPrInc,QjPrim,lPrim,lPrInc,QlPrim
            Write (6,*) MemMax,Mem0,Mem1,Mem2,Mem3
            Write (6,*) MemTrn,MemRys,MemPrm
            Call Abend()
         End If
         Go To 999
      End If
*-----Subtracte one additional word for getmem's internal error check
      Mem0 = Mem0 - Mem3 - 1
      MinXtr = Min(MinXtr,Mem0)
*
      MemSum=Mem1+Mem2+Mem3
      Mem2 = Mem2 + Mem3
*
      ipMem2 = ipMem1 + Mem1
*
      r1 = r1 + DBLE(iBsInc)/DBLE(iBas)
      r2 = r2 + DBLE(jBsInc)/DBLE(jBas)
      r3 = r3 + DBLE(kBsInc)/DBLE(kBas)
      r4 = r4 + DBLE(lBsInc)/DBLE(lBas)
      q1 = q1 + DBLE(iPrInc)/DBLE(iPrim)
      q2 = q2 + DBLE(jPrInc)/DBLE(jPrim)
      q3 = q3 + DBLE(kPrInc)/DBLE(kPrim)
      q4 = q4 + DBLE(lPrInc)/DBLE(lPrim)
*     Call GetMem('PSOAO1','CHECK','REAL',iDum,iDum)
*     Call qExit('PSOAO1')
      Return
      End
