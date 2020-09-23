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
* Copyright (C) 1990,2015, Roland Lindh                                *
*               1990, IBM                                              *
************************************************************************
      SubRoutine PSOAO0(nSO,MemPrm,MemMax,iAnga, iCmpa,
     &                  iBas,  iBsInc, jBas,  jBsInc,
     &                  kBas,  kBsInc, lBas,  lBsInc,
     &                  iPrim, iPrInc, jPrim, jPrInc,
     &                  kPrim, kPrInc, lPrim, lPrInc,
     &                  ipMem1,ipMem2,
     &                  Mem1,  Mem2,  DoFock)
************************************************************************
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
* Called from: Eval_Ints                                               *
*                                                                      *
* Calling    : QEnter                                                  *
*              Change                                                  *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*                                                                      *
*             Modified for scratch to FckAcc                           *
*                                                                      *
*             Modified for unified Work2 and Work3 block. Febr. 2015   *
************************************************************************
      use lw_Info
      use Temporary_parameters, only: force_part_c, force_part_p
      use Integral_Parameters, only: iWROpt
      use RICD_Info, only: Do_RI, Cholesky
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
#include "WrkSpc.fh"
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "lCache.fh"
#include "pstat.fh"
      Integer iAnga(4), iCmpa(4)
      Logical QiBas, QjBas, QkBas, QlBas, QjPrim, QlPrim, DoFock, Fail
#include "SysDef.fh"
*
*     Statement function to compute canonical index
*
      nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6  - 1
*
      la = iAnga(1)
      lb = iAnga(2)
      lc = iAnga(3)
      ld = iAnga(4)
      iCmp = iCmpa(1)
      jCmp = iCmpa(2)
      nab=iCmp*jCmp
      kCmp = iCmpa(3)
      lCmp = iCmpa(4)
      ncd=kCmp*lCmp
      iTotal = iTotal + 1
      mabMin=nabSz(Max(la,lb)-1)+1
      mabMax=nabSz(la+lb)
      ne=(mabMax-mabMin+1)
      mcdMin=nabSz(Max(lc,ld)-1)+1
      mcdMax=nabSz(lc+ld)
      nf=(mcdMax-mcdMin+1)
*     e0|f0 size
      mabcd=ne*nf
*     ab|cd size
      nabcd=nab*ncd
*
*     For test purposed we can force partitioning for both the
*     contracted and primitive block
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
      iFact = 1
      If (iWropt.eq.0 .and.
     &   .Not.(Cholesky.or.Do_RI)) iFact = 4 + 3
*
 999  Continue
*
*     We can partion the blocks on all contrcated indicies
*     and on the second and fourth primitive indicies.
*     We put priority in keeping an as large fraction of the contracted
*     block as much preserved to limit the recompution of primitive
*     integrals.
*
      QjPrim = .False.
      QlPrim = .True.
      QiBas  = .False.
      QjBas  = .False.
      QkBas  = .False.
      QlBas  = .False.
      Mem0 = MemMax
*
      ijPrInc=iPrInc*jPrInc
      klPrInc=kPrInc*lPrInc
      mijkl = iPrInc*jPrInc*kPrInc*lPrInc
      nijkl = iBsInc*jBsInc*kBsInc*lBsInc
*
*-----Work1
*     Memory for SO block. If petite list format is used there
*     will be no SO block.
*
      kSOInt = nSO*nijkl
      Mem1 = iFact*kSOInt
      If (Mem1.eq.0) Mem1 = 1
      If (nIrrep.eq.1) Mem1 = 1 + (iFact-1) * nabcd*nijkl
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
            Call WarningMessage(2,' Allocation failed for Work1')
            Write (6,*) Mem0,Mem1
            Write (6,*) iPrInc,iBsInc,kPrInc,kBsInc,
     &                  jPrInc,jBsInc,lPrInc,lBsInc
            Call Abend()
         End If
         Go To 999
      End If
      Mem0 = Mem0 - Mem1 - 1
*
*-----Work2
*     MemPr  : Scratch for Rys
*
*     Memory for the Rys-Gauss procedure. This includes memory for
*     all the intermediate arrays AND the {e0|f0} block.
*
      MemPr = MemPrm*mijkl
*
*     MemAux : Auxiliary memory for partial contraction.
*
*     If the primitive block is not full we need to accumulate primitive
*     contributions to the contracted block. This will require a
*     permanant pice of memory, during the computation of the primitive
*     sublocks, where the incomplete block of contracted integrals are
*     stored (Work4) - mabcd*nijkl.
*
      If (jPrInc.ne.jPrim.or.lPrInc.ne.lPrim) Then
         MemAux = Max(mabcd,nabcd)*nijkl
      Else
         MemAux = 0
      End If
*     Write (*,*) 'MemAux:',MemAux
*
*     MemCon : Scratch for the contraction step
*
*     We need memory to the largest of the {e0|f0} or (e0|f0) block.
*     Note that this correspond to that the A3 block in the second
*     halftransformation is set to zero. We also assume that the
*     integrals might be transformed to real spherical harmonic before
*     the contraction step. This only happens for mijkl>=nijl
*
*     The routine which does the contraction does this in two
*     half contraction steps with an intermediate which contains
*     a subset of the quater or three quarter transformed integrals.
*     Since we do not know the order of contraction at this time we will
*     have to compute for the worst case scenario.
*
*     Contraction of the two first indicies: iPrim->iBas & jPrim->jBas
*     while the third and fourth indicies are uncontrcated.
      nCache_ = (3*lCache)/4 - iPrim*iBas - jPrim*jBas
*     Note that we do not know here the order of the contraction, hence
*     we assume worst case scenario. The worst case is when lSize is
*     the smallest.
      lSize = iPrInc*jPrInc + Min(jPrInc*iBsInc,iPrInc*jBsInc)
*     Length of the compound index.
      nVec1 = kPrInc*lPrInc * Max(nabcd,mabcd)
*     Compute the subsection length of the compound index.
      IncVec = Min(Max(1,nCache_/lSize),nVec1)
*     Size of the initial integrals
      nA1a= Max(nabcd,mabcd)*mijkl
*     Size of the intemediate array
      nA2a= Max(iPrim,jPrim)*IncVec
*     Size of the half transformed integrals
      nA3a= iBsInc*jBsInc*nVec1
*     Write (*,*)
*     Write (*,*) 'IncVec,iPrim,jPrim:',IncVec,iPrim,jPrim
*     Write (*,*) 'nVec1,lSize=',nVec1,lSize
*     Write (*,*) 'nA1,nA2,nA3:',nA1a,nA2a,nA3a
*
*     Contraction of the two last indicies: kPrim->kBas & lPrim->lBas
*     while the first and second indicies are contrcated.
      nCache_ = (3*lCache)/4 - kPrim*kBas - lPrim*lBas
      lSize = kPrInc*lPrInc + Min(lPrInc*kBsInc,kPrInc*lBsInc)
      nVec2 = iBsInc*jBsInc * Max(nabcd,mabcd)
      IncVec = Min(Max(1,nCache_/lSize),nVec2)
      nA1b= nVec2*kPrInc*lPrInc
      nA2b= Max(kPrim,lPrim)*IncVec
      nA3b= kBsInc*lBsInc*nVec2
      If (MemAux.ne.0) nA3b=0
      MemCon=Max(nA1a,nA3b)+Max(nA2a,nA2b)+Max(nA3a,nA1b)
*     Write (*,*) 'IncVec,kPrim,lPrim:',IncVec,kPrim,lPrim
*     Write (*,*) 'nVec2,lSize=',nVec1,lSize
*     Write (*,*) 'nA1,nA2,nA3:',nA1b,nA2b,nA3b
*     Write (*,*) 'MemCon     :',MemCon
*
*     Contraction of the two last indicies: kPrim->kBas & lPrim->lBas
*     while the first and second indicies are uncontrcated.
      nCache_ = (3*lCache)/4 - kPrim*kBas - lPrim*lBas
      lSize = kPrInc*lPrInc + Min(lPrInc*kBsInc,kPrInc*lBsInc)
      nVec1 = iPrInc*jPrInc * Max(nabcd,mabcd)
      IncVec = Min(Max(1,nCache_/lSize),nVec1)
      nA1a= Max(nabcd,mabcd)*mijkl
      nA2a= IncVec*Max(kPrim,lPrim)
      nA3a= kBsInc*lBsInc*nVec1
*     Write (*,*) 'IncVec,kPrim,lPrim:',IncVec,kPrim,lPrim
*     Write (*,*) 'nVec1,lSize=',nVec1,lSize
*     Write (*,*) 'nA1,nA2,nA3:',nA1a,nA2a,nA3a
*
*     Contraction of the two first indicies: iPrim->iBas & jPrim->jBas
*     while the third and fourth indicies are contrcated.
      nCache_ = (3*lCache)/4 - iPrim*iBas - jPrim*jBas
      lSize = iPrInc*jPrInc + Min(jPrInc*iBsInc,iPrInc*jBsInc)
      nVec2 = kBsInc*lBsInc * Max(nabcd,mabcd)
      IncVec = Min(Max(1,nCache_/lSize),nVec2)
      nA1b= nVec2*iPrInc*jPrInc
      nA2b= IncVec*Max(iPrim,jPrim)
      nA3b= iBsInc*jBsInc*nVec2
      If (MemAux.ne.0) nA3b=0
      MemCon = Max(MemCon,Max(nA1a,nA3b)+Max(nA2a,nA2b)+Max(nA3a,nA1b))
*     Write (*,*) 'IncVec,iPrim,jPrim:',IncVec,iPrim,jPrim
*     Write (*,*) 'nVec2,lSize=',nVec1,lSize
*     Write (*,*) 'nA1,nA2,nA3:',nA1b,nA2b,nA3b
*     Write (*,*) 'MemCon     :',MemCon
*
*     MemSp1 : Scratch for the transformation from cartesian to
*              spherical harmonics. This is executed inside the
*              primitive block only if (1) the set of primitive
*              integrals are fewer than the contracted functions
*              and (2) there is no partitioning on the primitive
*              index.
*
      If (jPrInc.ne.jPrim.or.lPrInc.ne.lPrim) Then
         MemSp1= Max(mabcd+nab*nf,nabcd+nab*nf)*nijkl
      Else
         MemSp1= Max(mabcd+nab*nf,nabcd+nab*nf)*nijkl
      End If
*
*     MemFck : Scratch for FckAcc.
*
*     Memory to manupulate the 1-particle densities in FckAcc.
*
      If (DoFock) Then
         MemFck= 6*Max(iBsInc*lBsInc,iBsInc*kBsInc,jBsInc*lBsInc,
     &          jBsInc*kBsInc,iBsInc*jBsInc,kBsInc*lBsInc)
     &         + nijkl*nabcd
      Else
         MemFck=0
      End If
*
*     MemPck : Scratch for integral packing
*
*     Memory for integral packing.
*
      If (iWROpt.eq.0 .and. .Not.(Cholesky.or.Do_RI)) Then
         MemPck= 2 * nabcd*nijkl
      Else
         MemPck = 0
      End If
      Mem2 = Max((MemPr+MemAux),
     &           (MemCon+MemAux),
     &           (MemSp1+MemAux),
     &           MemFck,
     &           MemPck)
      If (Mem2+1.gt.Mem0) Then
         MaxReq=Max(MaxReq,Mem2+1-Mem0)
         Call Change(iBas, iBsInc,QiBas, kBas, kBsInc,QkBas,
     &               jBas, jBsInc,QjBas, lBas, lBsInc,QlBas,
     &               jPrim,jPrInc,QjPrim,lPrim,lPrInc,QlPrim,MaxReq,
     &               Fail)
         If (Fail) Then
            Call WarningMessage(2,' Allocation failed for Work2')
            Write (6,*) Mem0
            Write (6,*) iPrInc,iBsInc,kPrInc,kBsInc,
     &                  jPrInc,jBsInc,lPrInc,lBsInc
            Call Abend()
         End If
         Go To 999
      End If
      Mem0 = Mem0 - Mem2 - 1
*
      ipMem2 = ipMem1 + Mem1
*
      MemSum=Mem1+Mem2
*
*     Auxiliary memory for integral packing
*
      If (iWropt.eq.0 .and. .Not.(Cholesky.or.Do_RI)) Then
         If (nIrrep.eq.1) Then
            lPack = nabcd*nijkl
            lwInt = ipMem1
         Else
            lPack = kSOInt
            lwInt = ipMem1 + lPack
         End If
         lwSyB = lwInt + 2*lPack
         lwSqN = lwSyB + 2*lPack
      Else
         lPack = 0
         lwInt  = 0
         lwSyB  = 0
         lwSqN  = 0
      End If
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
      End
