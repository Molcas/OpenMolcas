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
* Copyright (C) 1990, Roland Lindh                                     *
*               1990, IBM                                              *
************************************************************************
      SubRoutine PSOAO0_h(nSO,nMemab,nMemcd,MemPrm,MemMax,
     &                            iAnga, iCmpa,
     &                            iBas,  iBsInc, jBas,  jBsInc,
     &                            kBas,  kBsInc, lBas,  lBsInc,
     &                            iPrim, iPrInc, jPrim, jPrInc,
     &                            kPrim, kPrInc, lPrim, lPrInc,
     &                            ipMem1,ipMem2,ipMem3,ipMem4, ipMend,
     &                              Mem1,  Mem2,  Mem3,  Mem4,   Mend)
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
* Called from: TwoEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              Change                                                  *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
************************************************************************
      use Integral_Parameters, only: iWROpt
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
#include "WrkSpc.fh"
#include "real.fh"
#include "print.fh"
#include "lCache.fh"
#include "pstat.fh"
#include "warnings.fh"
      Integer iAnga(4), iCmpa(4)
      Logical QiBas, QjBas, QkBas, QlBas, QjPrim, QlPrim, Fail
#include "SysDef.fh"
*
*     Statement function to compute canonical index
*
      nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6  - 1
*
      iRout = 10
      iPrint = nPrint(iRout)
*     iQ = 1
*     Call qEnter('PSOAO0')
      la = iAnga(1)
      lb = iAnga(2)
      lc = iAnga(3)
      ld = iAnga(4)
      iCmp = iCmpa(1)
      jCmp = iCmpa(2)
      kCmp = iCmpa(3)
      lCmp = iCmpa(4)
      iTotal = iTotal + 1
      mabMin=nabSz(Max(la,lb)-1)+1
      mabMax=nabSz(la+lb)
      mcdMin=nabSz(Max(lc,ld)-1)+1
      mcdMax=nabSz(lc+ld)
      mabcd=(mabMax-mabMin+1)*(mcdMax-mcdMin+1)
      nabcd=iCmp*jCmp*kCmp*lCmp
*
      iBsInc = iBas
      jBsInc = jBas
      kBsInc = kBas
      lBsInc = lBas
      iPrInc = iPrim
      jPrInc = jPrim
      kPrInc = kPrim
      lPrInc = lPrim
      iFact = 1
      If (iWropt.eq.0) iFact = 4/RtoI+3
*
 999  Continue
      QjPrim = .False.
      QlPrim = .True.
      QiBas  = .False.
      QjBas  = .False.
      QkBas  = .False.
      QlBas  = .False.
      Mem0 = MemMax
*
*-----Work1
*     Memory for SO block. If petite list format is used there
*     will be no SO block.
*
      kSOInt = nSO*iBsInc*jBsInc*kBsInc*lBsInc
      Mem1 = iFact*kSOInt
      If (Mem1.eq.0) Mem1 = 1
      If (nIrrep==1) Mem1 = 1 + (iFact-1) *
     &                   iCmp*  jCmp*  kCmp*  lCmp*
     &                   iBsInc*jBsInc*kBsInc*lBsInc
      If (Mem1+1.gt.Mem0) Then
         MaxReq=Max(MaxReq,Mem1+1-Mem0)
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
            Write (6,*) ' Allocation failed for Work1'
            Write (6,*) Mem0,Mem1
            Write (6,*) iPrInc,iBsInc,kPrInc,kBsInc,
     &                  jPrInc,jBsInc,lPrInc,lBsInc
            Call Quit(_RC_MEMORY_ERROR_)
         End If
         Go To 999
      End If
*     Write (*,*) ' After Mem1', iPrInc,iBsInc,kPrInc,kBsInc,
*    &                           jPrInc,jBsInc,lPrInc,lBsInc
      Mem0 = Mem0 - Mem1 - 1
*
*-----Work2
*     MemPr  : Scratch for Rys
*     MemCon : Scratch for the contraction step
*     MemTr1 : Scratch for the 1st application of the transfer eqn.
*     MemTr2 : Scratch for the 2nd application of the transfer eqn.
*-----Work4 (this is overlayed with Work2 and is placed at the top)
*     MemAux : Auxiliary memory for partial contraction.
*-----Work5 (this is overlayed with Work2 and is placed at the top)
*
      MemPr = MemPrm*iPrInc*jPrInc*kPrInc*lPrInc
      MemCon= mabcd*Max(iPrInc*jPrInc*kPrInc*lPrInc,
     &                  iBsInc*jBsInc*kBsInc*lBsInc)
      If (jPrInc.ne.jPrim.or.lPrInc.ne.lPrim) Then
         MemAux = mabcd*iBsInc*jBsInc*kBsInc*lBsInc
      Else
         MemAux = 0
      End If
      MemTr1= (mabMax-mabMin+1)*nMemcd*iBsInc*jBsInc*kBsInc*lBsInc
      MemTr2= kCmp*lCmp*nMemab*iBsInc*jBsInc*kBsInc*lBsInc
      Mem2  = Max(MemPr+MemAux,MemCon+MemAux,MemTr1,MemTr2)
      If (Mem2+1.gt.Mem0) Then
         MaxReq=Max(MaxReq,Mem2+1-Mem0)
         Call Change(iBas, iBsInc,QiBas, kBas, kBsInc,QkBas,
     &               jBas, jBsInc,QjBas, lBas, lBsInc,QlBas,
     &               jPrim,jPrInc,QjPrim,lPrim,lPrInc,QlPrim,MaxReq,
     &               Fail)
         If (Fail) Then
            Write (6,*) ' Allocation failed for Work2'
            Write (6,*) Mem0,Mem2,MemPr+MemAux,MemCon+MemAux,
     &                  MemTr1,MemTr2
            Write (6,*) iPrInc,iBsInc,kPrInc,kBsInc,
     &                  jPrInc,jBsInc,lPrInc,lBsInc
            Call Quit(_RC_MEMORY_ERROR_)
         End If
         Go To 999
      End If
      If (jPrInc.ne.jPrim.or.lPrInc.ne.lPrim) Then
         Mem4 = MemAux
      Else
         Mem4 = Mem2
      End If
*     Write (*,*) ' After Mem2', jPrInc,jBsInc,lPrInc,lBsInc
      Mem0 = Mem0 - Mem2 - 1
*
*-----Work3
*     MemCon : Scratch for the contraction step, and transpose after the
*              contraction step
*     MemSp1 : Scratch for the 1st transformation from cartesian to
*              spherical harmonics.
*     MemSp2 : Scratch for the 2nd transformation from cartesian to
*              spherical harmonics.
*     MemTr3 : Scratch for transpose in tnsctl
*
      nCache_ = (3*lCache)/4 - iPrim*iBas - jPrim*jBas
      lSize = iPrInc*jPrInc + Max(jPrInc*iBsInc,iPrInc*jBsInc)
      nVec1 = kPrInc*lPrInc * mabcd
      IncVec = Min(Max(1,nCache_/lSize),nVec1)
      nA2 = Max(jPrInc*iBsInc,iPrInc*jBsInc)*IncVec
      nA3 = iBsInc*jBsInc*nVec1
      MemCon = Max(mabcd*iBsInc*jBsInc*kBsInc*lBsInc,nA2+nA3)
*
      nCache_ = (3*lCache)/4 - kPrim*kBas - lPrim*lBas
      lSize = kPrInc*lPrInc + Max(lPrInc*kBsInc,kPrInc*lBsInc)
      nVec2 = iBsInc*jBsInc * mabcd
      IncVec = Min(Max(1,nCache_/lSize),nVec2)
      nA2 = Max(lPrInc*kBsInc,kPrInc*lBsInc)*IncVec
*     nA3 = kBsInc*lBsInc*nVec2
      MemCon = Max(MemCon,nA3+nA2)
*
************************************************************************
*
      nCache_ = (3*lCache)/4 - kPrim*kBas - lPrim*lBas
      lSize = kPrInc*lPrInc + Max(lPrInc*kBsInc,kPrInc*lBsInc)
      nVec1 = iPrInc*jPrInc * mabcd
      IncVec = Min(Max(1,nCache_/lSize),nVec1)
      nA2 = IncVec*Max(lPrInc*kBsInc,kPrInc*lBsInc)
      nA3 = kBsInc*lBsInc*nVec1
      MemCon = Max(MemCon,nA3+nA2)
*
      nCache_ = (3*lCache)/4 - iPrim*iBas - jPrim*jBas
      lSize = iPrInc*jPrInc + Max(jPrInc*iBsInc,iPrInc*jBsInc)
      nVec2 = kBsInc*lBsInc * mabcd
      IncVec = Min(Max(1,nCache_/lSize),nVec2)
      nA2 = IncVec*Max(jPrInc*iBsInc,iPrInc*jBsInc)
*     nA3 = iBsInc*jBsInc*nVec2
      MemCon = Max(MemCon,nA3+nA2)
*
      MemSp1 = (mabMax-mabMin+1)*lCmp*(lc+1)*(lc+2)/2 *
     &       iBsInc*jBsInc*kBsInc*lBsInc
      MemSp2= lCmp*kCmp*jCmp*(la+1)*(la+2)/2 *
     &       iBsInc*jBsInc*kBsInc*lBsInc
      MemTr3 = mabcd*iBsInc*jBsInc*kBsInc*lBsInc
      Mem3 = Max(MemCon,MemSp1,MemSp2,MemTr3)
      If (Mem3+1.gt.Mem0) Then
         MaxReq=Max(MaxReq,Mem3+1-Mem0)
         Call Change(iBas, iBsInc,QiBas, kBas, kBsInc,QkBas,
     &               jBas, jBsInc,QjBas, lBas, lBsInc,QlBas,
     &               jPrim,jPrInc,QjPrim,lPrim,lPrInc,QlPrim,MaxReq,
     &               Fail)
         If (Fail) Then
            Write (6,*) ' Allocation failed for Work3'
            Write (6,*) Mem0,Mem3,MemCon,MemSp1,MemSp2
            Write (6,*) iPrInc,iBsInc,kPrInc,kBsInc,
     &                  jPrInc,jBsInc,lPrInc,lBsInc
            Call Quit(_RC_MEMORY_ERROR_)
         End If
         Go To 999
      End If
      Mem0 = Mem0 - Mem3 - 1
      MinXtr = Min(MinXtr,Mem0)
*
      ipMem2 = ipMem1 + Mem1
      ipMem3 = ipMem2 + Mem2
      ipMem4 = ipMem2 + Mem2 - Mem4
      Mend=0
*
      MemSum=Mem1+Mem2+Mem3
*
      r1 = r1 + DBLE(iBsInc)/DBLE(iBas)
      r2 = r2 + DBLE(jBsInc)/DBLE(jBas)
      r3 = r3 + DBLE(kBsInc)/DBLE(kBas)
      r4 = r4 + DBLE(lBsInc)/DBLE(lBas)
      q1 = q1 + DBLE(iPrInc)/DBLE(iPrim)
      q2 = q2 + DBLE(jPrInc)/DBLE(jPrim)
      q3 = q3 + DBLE(kPrInc)/DBLE(kPrim)
      q4 = q4 + DBLE(lPrInc)/DBLE(lPrim)
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(ipMend)
      End
