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
*               1995, Anders Bernhardsson                              *
************************************************************************
      SubRoutine PSOAO2(nSO,MemPrm,MemM,
     &                            iAnga, iCmpa, iShela, iFnc,
     &                            iBas,  iBsInc, jBas,  jBsInc,
     &                            kBas,  kBsInc, lBas,  lBsInc,
     &                            iPrim, iPrInc, jPrim, jPrInc,
     &                            kPrim, kPrInc, lPrim, lPrInc,
     &                            nAco,
     &                            Mem1,Mem2,Mem3,Mem4,
     &                            MemX,MemPSO,
     &                            MemFck,nFT,nCMO,
     &                            MemFin,MemBuffer,
     &                            iMemB)
************************************************************************
*                                                                      *
*  Object: to partion the SO and AO block. It will go to some length   *
*          before it will start and break up the SO block. This will   *
*          reduce the total flop count. However, as we are breaking up *
*          the AO block this will affect the vectorization. Hence, at  *
*          some point it will actually be better to recompute the      *
*          primitives.                                                 *
*          Current stratergy:                                          *
*          1. Reduce the size of the density matrix and buffer so that *
*             it fits into memory.                                     *
*                                                                      *
*          2. Start reducing the length of the primitives in the order *
*             lPrim,jPrim.                                             *
*                                                                      *
*          3. Reduce the size of the SO block by reducing the number of*
*             basis functions in the order lBas, jBas.                 *
*                                                                      *
*          4. Reduce the size of the Buffer.                           *
*                                                                      *
*          5. Reduce kBas,iBas                                         *
*                                                                      *
*          6. Terminate run telling job max and min of additional      *
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
*             Modified to first order derivatives. January '92         *
*             Anders Bernhardsson Theoretical chemistry, Lund 1995     *
************************************************************************
*
*               Memory map for mckinley
*
*---------------------------------------------------------------------------
*|      |               |       |               |               |          |
*|REAL  |  P TRANSF     | RYSG2 |   Transf      | FCK GENERAT   |MO Transf |
*|      |               |       |               |               |          |
*---------------------------------------------------------------------------
*|      |               |       |               |               |Scratch   |
*|  MX  |               |       | 9*abcd*ijkl   |    SS         |space     |
*|      |               |       |               |               |          |
*---------------------------------------------------------------------------
*|      |               |       |               |               |          |
*|  M3  |Scratch space  |Memrys |Scratch space  |    SS         | Scratch  |
*|      |               |       |               |               | space    |
*---------------------------------------------------------------------------
*|      |MEM4 (half tr) |*******|***************|               |          |
*|  M2  |               |       |               |               |  SS      |
*|      |PSO transf     |       |               |Scratch space  |          |
*---------------------------------------------------------------------------
*|      |               |       |               |               |          |
*|  M1  |      P        |   *   |      *        |     *         |    *     |
*|      |               |       |               |               |          |
*---------------------------------------------------------------------------
*|      |      ?        |       |               |               |          |
*|Buffer|***************|*******|Transformed    |***************|**********|
*|      |               |       |integrals      |               |          |
*---------------------------------------------------------------------------
*
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
c#include "print.fh"
#include "pstat.fh"
#include "pso.fh"
#include "disp.fh"
#include "disp2.fh"
#include "buffer.fh"
      Integer iAnga(4), iCmpa(4), nPam(4,0:7), iiBas(4), iShela(4),
     &        iFnc(4)
      Logical QiBas, QjBas, QkBas, QlBas, QjPrim, QlPrim, Fail
      Integer iTwoj(0:7),iMemB
      Data iTwoj/1,2,4,8,16,32,64,128/
*
*     Statement function to compute canonical index
*
      nElem(i) = (i+1)*(i+2)/2
*
c     iRout = 10
c     iPrint = nPrint(iRout)
c     Call qEnter('PSOAO2')
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
*     If (force_part_c) Then
*        iBsInc = (iBas+1)/2
*        jBsInc = (jBas+1)/2
*        kBsInc = (kBas+1)/2
*        lBsInc = (lBas+1)/2
*     Else
         iBsInc = iBas
         jBsInc = jBas
         kBsInc = kBas
         lBsInc = lBas
*     End If
      If (force_part_p) Then
         jPrInc = (jPrim+1)/2
*        lPrInc = (lPrim+1)/2
      Else
         jPrInc = jPrim
*        lPrInc = lPrim
      End If
      iPrInc = iPrim
      kPrInc = kPrim
      lPrInc = lPrim
      MemBuffer = iMemB
      MemMax=MemM-MemBuffer
*
 999  Continue
      nijkl=iBsInc*jBsInc*kBsInc*lBsInc
      QjPrim = .False.
      QlPrim = .True.
      QiBas  = .False.
      QjBas  = .False.
      QkBas  = .False.
      QlBas  = .False.
      Mem0 = MemMax
*
*    Picked MO coeff
*
      If (nMethod.eq.RASSCF) Then
       nCMO=nACO*kCmp*kBas+nACO*lCmp*lBas
      Else
       nCMO=0
      End If
*
*     Area for integral storage before transforming them to FM/MO
*     and place for the CMOs
*
      MemFin=9*nijkl*nabcd
      If (MemFin+ncmo+1.gt.Mem0) Then
         MaxReq=Max(MaxReq,nCMO+MemFin+1-Mem0)
         QlPrim=.false.
         Call Change(iBas, iBsInc,QiBas, kBas, kBsInc,QkBas,
     &               jBas, jBsInc,QjBas, lBas, lBsInc,QlBas,
     &               jPrim,jPrInc,QjPrim,lPrim,lPrInc,QlPrim,MaxReq,
     &               Fail)
         If (Fail) Then
            Write (6,*) 'PSOAO2: memory partitioning failed!'
            Write (6,*) '        Restart with more memory!'
            Call QTrace
            Call Abend()
         End If
         Go To 999
      End If
*-----Subtract one additional word for getmem's internal error check
      Mem0 = Mem0 - MemFin -nCMO - 1
*
*
*--------------------------------------------------------------------
*
*-----*** Work1 ***
*
*
*-----Memory for 2nd order density matrix in SO basis.
*
      kSOInt = nSO*iBsInc*jBsInc*kBsInc*lBsInc
      Mem1 = kSOInt
*
*-----Allocate memory for MO to SO/AO transformation
*     of the 2nd order density matrix for this shell quadruplet.
*     and area for AO/SO transformation of Fock matrix.
*
      If (lPSO) Then
         iiBas(1) = iBsInc
         iiBas(2) = jBsInc
         iiBas(3) = kBsInc
         iiBas(4) = lBsInc
         Call ICopy(4*8,[0],0,nPam,1)
         MemPSO = 1
         nTmp2 = 0
*        Call IecPrt('iiBas',iiBas,1,4)
*
         Do 9 jPam = 1, 4
            iTmp1= 0
            nTmp1= 0
            Do 10 j = 0, nIrrep-1
               Do 11 i1 = 1, iCmpa(jPam)
                  If (iAnd(IrrCmp(IndS(iShela(jPam))+i1),
     &                iTwoj(j)).ne.0) Then
                      nPam(jPam,j) = nPam(jPam,j) + iiBas(jPam)
                      nTmp1= nTmp1+ iiBas(jPam)
                      iTmp1= iTmp1+ 1
                  End If
 11            Continue
 10         Continue
            MemPSO = MemPSO * nTmp1
            nTmp2 = nTmp2 + nTmp1
            iFnc(jPam) = iTmp1
 9       Continue
         MemScr=MemTra(nPam)
         nFac = 4
         nTmp2 = nTmp2 + 4
      Else
         MemScr=0
         MemPSO=0
         nFac = 0
         nTmp2 = 0
      End If
      MemAux = MemPSO + MemScr + nFac*nDim + nTmp2 + 4
      MemB_AUX=MemAux
      If (Mem1+1+MemAux.gt.Mem0) Then
         MaxReq=Max(MaxReq,Mem1+1+MemAux-Mem0)
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
            Write (6,*) 'PSOAO2: memory partitioning failed!'
            Write (6,*) '        Restart with more memory!'
            Call QTrace
            Call Abend()
         End If
         Go To 999
      End If
      Mem0 = Mem0 - Mem1 - 1
*---------------------------------------------------------------
*
*     MemFck: Target for generating the symmetrized Fock Matrix.
*              Distributed localy.
*     MemMo : Target for generating the  MO integrals
*              Distributed localy.
*    Whole work area is used, if work area is big enough
*    work3 is increased
*
*---------------------------------------------------------------
*
      MemDep=nijkl*nabcd
*     Temp+S1+S2
      MemFck=2*MemDep+Max(MemDep,nijkl+
     &        Max(iBsInc*lBsInc,jBsInc*lBsInc,
     &            iBsInc*kBsInc,jBsInc*kBsInc))
      nFT=iBsInc*jBsInc*iCmp*jCmp+kBsInc*lBsInc*kCmp*lCmp+
     &    iBsInc*kBsInc*iCmp*kCmp+jBsInc*lBsInc*jCmp*lCmp+
     &    iBsInc*lBsInc*iCmp*lCmp+jBsInc*kBsInc*jCmp*kCmp
      MemFck=MemFck+nFT
      If (nmethod.eq.RASSCF) Then
*
*     3 scratch spaces, sorted integrals and translation
*
       nMaxC=nACO
       MemFck=MemFck+2*nMaxC
       nMax=Max(iCmp*iBsInc,jCmp*jBsInc,kCmp*kBsInc,lcmp*lBsInc)
       nMax=Max(nMax,nMaxC)
       memMO=3*nMax**4+10*nabcd*nijkl
      Else
       MemMo=0
      End If
*
*---------------------------------------------------------------------
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

      MemDeP = nabcd * nijkl
      MemTrn = mabcd * Max(iBsInc*jBsInc*kBsInc*lBsInc,
     &                     iPrInc*jPrInc*kBsInc*lBsInc,
     &                     iPrInc*jPrInc*kPrInc*lPrInc)
      MemTrn=MemTrn+1

*-----If partial decontraction we need to keep the contracted 2nd
*     order density matrix. (Work4)
      If (jPrInc.ne.jPrim .or. lPrInc.ne.lPrim) Then
         MemAux = mabcd*iBsInc*jBsInc*kBsInc*lBsInc
      Else
         MemAux = 0
      End If
      MemSph = mabcd * iBsInc*jBsInc*kBsInc*lBsInc
      Mem2 = Max(MemTrn+MemAux,MemDeP,MemSph)
      MemFck=MemFck-Mem2
      MemMO=MemMo-Mem2
      If (Mem2+1.gt.Mem0) Then
         MaxReq=Max(MaxReq,Mem2+1-Mem0)
         Call Change(iBas, iBsInc,QiBas, kBas, kBsInc,QkBas,
     &               jBas, jBsInc,QjBas, lBas, lBsInc,QlBas,
     &               jPrim,jPrInc,QjPrim,lPrim,lPrInc,QlPrim,MaxReq,
     &               Fail)
         If (Fail) Then
            Write (6,*) 'PSOAO2: memory partitioning failed!'
            Write (6,*) '        Restart with more memory!'
            Call QTrace
            Call Abend()
         End If
         Go To 999
      End If
*-----Subtracte one additional word for getmem's internal error check
      Mem0 = Mem0 - Mem2 - 1
      MemX = 9*mabcd * iBsInc*jBsInc*kBsInc*lBsInc
      MemFck=MemFck-MemX
      MemMO=MemMo-MemX
      If (MemX+1.gt.Mem0) Then
         MaxReq=Max(MaxReq,MemX+1-Mem0)
         Call Change(iBas, iBsInc,QiBas, kBas, kBsInc,QkBas,
     &               jBas, jBsInc,QjBas, lBas, lBsInc,QlBas,
     &               jPrim,jPrInc,QjPrim,lPrim,lPrInc,QlPrim,MaxReq,
     &               Fail)
         If (Fail) Then
            Write (6,*) 'PSOAO2: memory partitioning failed!'
            Write (6,*) '        Restart with more memory!'
            Call QTrace
            Call Abend()
         End If
         Go To 999
      End If
      Mem0=Mem0-MemX

*
*-----*** Work3 and Work5 ***
*
*-----Scratch for decontraction and transformation to spherical gaussian.
*     Working array for Rysg2.
*     Scratch area for resolving degeneracies due to the double coset
*     treatement of the symmetry.
*     MemTrn: Scratch for decontraction
*     MemRys: Scratch for calualation of primitive integral gradients.
*
      iFac = 1
      If (mabcd.ne.1) iFac = 2
      MemF=9*nabcd*nijkl
      MemTrn=mabcd * Max(iPrInc*jBsInc*kBsInc*lBsInc,
     &                   iPrInc*jPrInc*kPrInc*lBsInc,
     &                   iPrInc*jPrInc*kPrInc*lPrInc*iFac)
      MemRys=MemPrm * iPrInc*jPrInc*kPrInc*lPrInc+80
*
*  Scratch space for contraction of the integrals
*
      MemCntrct=9*mabcd*(Max(iBsInc*jBsInc*kBsInc*lPrInc,
     &                       iBsInc*jPrInc*kPrInc*lPrInc)+
     &                       iBsInc*jBsInc*kPrInc*lPrInc)

      MemFck=Max(0,MemFck)
      MemMo=Max(0,MemMo)
      Mem3 = Max(MemMO,MemFck,MemTrn, MemRys, 2*MemF,
     &           MemF+MemCntrct)
      If (Mem3+1.gt.Mem0) Then
         MaxReq=Max(MaxReq,Mem3+1-Mem0)
         Call Change(iBas, iBsInc,QiBas, kBas, kBsInc,QkBas,
     &               jBas, jBsInc,QjBas, lBas, lBsInc,QlBas,
     &               jPrim,jPrInc,QjPrim,lPrim,lPrInc,QlPrim,MaxReq,
     &               Fail)
         If (Fail) Then
            Write (6,*) 'PSOAO2: memory partitioning failed!'
            Write (6,*) '        Restart with more memory!'
            Call QTrace
            Call Abend()
         End If
         Go To 999
      End If
*-----Subtracte one additional word for getmem's internal error check
      Mem0 = Mem0 - Mem3 - 1
*
*-----Work4, if used, is placed at the end of Work2
      If (jPrInc.ne.jPrim .or. lPrInc.ne.lPrim) Then
         Mem4 = MemAux
      Else
         Mem4 = Mem2
      End If
*
      MemSum=Mem1+Mem2+Mem3+MemX+MemFin
c     Call qExit('PSOAO2')
      Return
      End
