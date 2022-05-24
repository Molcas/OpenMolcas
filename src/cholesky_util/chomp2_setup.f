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
* Copyright (C) 2004,2005, Thomas Bondo Pedersen                       *
************************************************************************
      SubRoutine ChoMP2_Setup(irc)
      use ChoMP2
C
C     Thomas Bondo Pedersen, Oct. 2004 / Feb. 2005.
C
C     Purpose: setup of Cholesky MP2 program.
C
#include "implicit.fh"
#include "cholesky.fh"
#include "choorb.fh"
#include "chomp2_cfg.fh"
#include "chomp2.fh"
#include "WrkSpc.fh"

      Character*12 SecNam
      Parameter (SecNam='ChoMP2_Setup')

      Character*2 Unt

      Integer NumVec(8), nFrac(8)
      Logical Accepted

      Integer bsize, blast
      Real*8  lX

      Logical ChoMP2_Setup_MemChk

      MulD2h(i,j)=iEor(i-1,j-1) + 1

      irc = 0

C     Setup index arrays and counters.
C     --------------------------------

      If ((DecoMP2 .or. DoDens) .and. ThrMP2.le.0.0D0) Then
         Call Get_dScalar('Cholesky Threshold',ThrMP2)
      End If

      Call ChoMP2_GetInf(nOrb,nOcc,nFro,nDel,nVir)
      iOcc(1) = 0
      iBatOrb(1) = 0
      iVir(1) = 0
      iFro(1) = 0
      iDel(1) = 0
      nOccT = nOcc(1)
      nVirT = nVir(1)
      nFroT = nFro(1)
      nDelT = nDel(1)
      If(.false.) Then
         nBatOrbT = nVir(1) + nOcc(1) + nFro(1) + nDel(1)
      Else
         nBatOrbT = nOcc(1)
      End If
      Do iSym = 2,nSym
         iOcc(iSym) = nOccT
         iVir(iSym) = nVirT
         iFro(iSym) = nFroT
         iDel(iSym) = nDelT
         iBatOrb(iSym) = nBatOrbT
         nOccT = nOccT + nOcc(iSym)
         nVirT = nVirT + nVir(iSym)
         nFroT = nFroT + nFro(iSym)
         nDelT = nDelT + nDel(iSym)
         If(.false.) Then
            nBatOrbT = nBatOrbT + nOcc(iSym) + nVir(iSym) + nFro(iSym)
     &                          + nDel(iSym)
         Else
            nBatOrbT = nBatOrbT + nOcc(iSym)
         End If
      End Do
*
      Do iSym = 1,nSym
         nT1am(iSym) = 0
         Do iSymi = 1,nSym
            iSyma = MulD2h(iSymi,iSym)
            iT1am(iSyma,iSymi) = nT1am(iSym)
            nT1am(iSym) = nT1am(iSym)
     &                  + nVir(iSyma)*nOcc(iSymi)
         End Do
      End Do
*
      If(.false.) Then
         Do iSym = 1, nSym
            nPQ_prod(iSym) = 0
            nPQ_prodij(iSym) = 0
            nPQ_prodia(iSym) = 0
            nPQ_prodab(iSym) = 0
            Do iSymQ = 1,nSym
               iSymP = MulD2h(iSymQ,iSym)
               iPQ_prod(iSymP,iSymQ) = nPQ_prod(iSym)

               nPQ_prod(iSym) = nPQ_Prod(iSym) +
     &                          (nOcc(iSymP)+nVir(iSymP)
     &                        +  nFro(iSymP)+nDel(iSymP))
     &                        * (nOcc(iSymQ)+nVir(iSymQ)
     &                        +  nFro(iSymQ)+nDel(iSymQ))
               nPQ_prodij(iSym) = nPQ_Prodij(iSym) +
     &                          (nFro(iSymP)+nOcc(iSymP))
     &                        * (nFro(iSymQ)+nOcc(iSymQ))
               nPQ_prodia(iSym) = nPQ_Prodia(iSym) +
     &                          (nFro(iSymP)+nOcc(iSymP))
     &                        * (nVir(iSymQ)+nDel(iSymQ))
               nPQ_prodab(iSym) = nPQ_Prodab(iSym) +
     &                          (nVir(iSymP)+nDel(iSymP))
     &                        * (nVir(iSymQ)+nDel(iSymQ))
            End Do
         End Do

      End If


      Do iSym = 1,nSym
         nT1AOT(iSym) = 0
         Do iSymAl = 1,nSym
            iSymi = MulD2h(iSymAl,iSym)
            iT1AOT(iSymi,iSymAl) = nT1AOT(iSym)
            nT1AOT(iSym) = nT1AOT(iSym)
     &                   + nOcc(iSymi)*nBas(iSymAl)
         End Do
      End Do

      Do iSym = 1,nSym
         nAOVir(iSym) = 0
         Do iSyma = 1,nSym
            iSymAl = MulD2h(iSyma,iSym)
            iAOVir(iSymAl,iSyma) = nAOVir(iSym)
            nAOVir(iSym) = nAOVir(iSym)
     &                   + nBas(iSymAl)*nVir(iSyma)
         End Do
      End Do

      If (ChoAlg .eq. 2) Then
         Do iSym = 1,nSym
            nMatab(iSym) = 0
            Do iSymb = 1,nSym
               iSyma = MulD2h(iSymb,iSym)
               iMatab(iSyma,iSymb) = nMatab(iSym)
               nMatab(iSym) = nMatab(iSym) + nVir(iSyma)*nVir(iSymb)
            End Do
         End Do
      Else
         Call Cho_iZero(nMatab,8)
         Call Cho_iZero(iMatab,64)
      End If

C     If batching over occuped orbitals is forced by user, there better
C     be orbitals to batch over!
C     -----------------------------------------------------------------

      If (ForceBatch .and. nBatOrbT.eq.1) ForceBatch = .false.

C     Setup batches over occupied orbitals.
C     -------------------------------------

*     nPQ_prodx will be the largest number of products in one symmetry
      nPQProdx = 0
      If(.false.) Then
         nPQProdx = nPQ_Prod(1)
         Do iSym = 2,nSym
            nPQProdx = max(nPQProdx,nPQ_Prod(iSym))
         End Do
      End If
*
      nT1amx = nT1am(1)
      Do iSym = 2,nSym
         nT1amx = max(nT1amx,nT1am(iSym))
      End Do

      If (DecoMP2) Then
         Do iSym = 1,nSym
            NumVec(iSym) = min(NumCho(iSym),nT1am(iSym))
         End Do
      Else
         Call iCopy(nSym,NumCho,1,NumVec,1)
      End If
      Call Cho_GAiGOp(NumVec,nSym,'max')

*     nBatOrbT is the total numbers of orbitals to batch over,
*     if densities are not to be computed this is set equal to nOccT

      If (nBatOrbT .lt. 6) Then
         mBatch = nBatOrbT
      Else
         mBatch = nBatOrbT/2 + 1
      End If
      Accepted = .false.
      nBatch = 0
      Do While (nBatch.lt.mBatch .and. .not.Accepted)
*
         nBatch = nBatch + 1
         If (nBatch .eq. mBatch) Then
            nBatch = nBatOrbT
            Do iSym = 1,nSym
               nFrac(iSym) = max(NumVec(iSym),1)
            End Do
         Else
            Do iSym = 1,nSym
               nFrac(iSym) = max(min(NumVec(iSym),20),1)
            End Do
         End If
*
         l_First     = nBatch
         l_FirstS    = nSym*nBatch
         l_NumOcc    = nBatch
         l_LnOcc     = nSym*nBatch
         l_LnT1am    = nSym*nBatch
         l_LiT1am    = nSym*nSym*nBatch
         l_NumBatOrb = nBatch
         l_LnBatOrb  = nSym*nBatch
         If (ChoAlg .eq. 2) Then
            l_LnMatij = nSym*nBatch
            l_LiMatij = nSym*nSym*nBatch
         Else
            l_LnMatij = 1
            l_LiMatij = 1
         End If
         If(.false.) Then
            l_LnPQprod = nSym*nBatch
            l_LiPQprod = nSym*nSym*nBatch
         Else
            l_LnPQprod = 1
            l_LiPQprod = 1
         End If
         l_lUnit  = nSym*nBatch
         Call ChoMP2_deallocate(irc)
         ChoMP2_allocated=.TRUE.
         Call GetMem('First','Allo','Inte',ip_First,l_First)
         Call GetMem('FirstS','Allo','Inte',ip_FirstS,l_FirstS)
         Call GetMem('NumOcc','Allo','Inte',ip_NumOcc,l_NumOcc)
         Call GetMem('LnOcc','Allo','Inte',ip_LnOcc,l_LnOcc)
         Call GetMem('LnT1am','Allo','Inte',ip_LnT1am,l_LnT1am)
         Call GetMem('LiT1am','Allo','Inte',ip_LiT1am,l_LiT1am)
         Call GetMem('LnMatij','Allo','Inte',ip_LnMatij,l_LnMatij)
         Call GetMem('LiMatij','Allo','Inte',ip_LiMatij,l_LiMatij)
         Call GetMem('lUnit','Allo','Inte',ip_lUnit,l_lUnit)
*     Generalization of NumOcc for arbitrary quantity to batch over
*     Would be good to kill NumOcc safely and only use one...
         Call GetMem('NumBatOrb','Allo','Inte',ip_NumBatOrb,l_NumBatOrb)
*     Generalization of LnOcc for arbitrary quantity to batch over.
         Call GetMem('LnBatOrb','Allo','Inte',ip_LnBatOrb,l_LnBatOrb)
         Call GetMem('LnPQprod','Allo','Inte',ip_LnPQprod,l_LnPQprod)
         Call GetMem('LiPQprod','Allo','Inte',ip_LiPQprod,l_LiPQprod)
         Call ChoMP2_Setup_Index(iWork(ip_First),iWork(ip_FirstS),
     &                           iWork(ip_NumOcc),iWork(ip_LnOcc),
     &                           iWork(ip_NumBatOrb),iWork(ip_LnBatOrb),
     &                           iWork(ip_LnT1am),iWork(ip_LiT1am),
     &                           iWork(ip_LnPQprod),iWork(ip_LiPQprod),
     &                           iWork(ip_LnMatij),iWork(ip_LiMatij),
     &                           nSym,nBatch)

         Call GetMem('MaxMP2','Max ','Real',ip_Dum,lWork)
         If(.false.) Then
*           All Memory available minus one full vector and some small
*           vectors for the PCG-algorithm.
            lAvail = lWork - nPQprodx-l_Mp2Lagr*9
         Else If (Laplace .and. SOS_MP2) Then
            lX=0.0d0
            Do iSym=1,nSym
               If (nT1am(iSym).gt.0 .and. NumVec(iSym).gt.0) Then
                  bsize=min(Laplace_BlockSize,NumVec(iSym))
                  nBlock=(NumVec(iSym)-1)/bsize+1
                  blast=NumVec(iSym)-bsize*(nBlock-1)
                  xM=dble(NumVec(iSym))
                  xn=dble(nBlock)
                  xb=dble(bsize)
                  xbp=dble(blast)
                  lX=max(lX,0.5d0*(xM*(xM+1.0d0)
     &                            +(xn-1.0d0)*xb*(xb-1.0d0)
     &                            +xbp*(xbp-1.0d0)))
               End If
            End Do
            l_X=int(lX)
            If (l_X .lt. 0) Then
               Write(Lupri,'(A,A)')
     &         SecNam,': dimension of X matrix is negative!'
               Write(Lupri,'(A,I15)') 'l_X=',l_X
               If (lX .gt. 0.0d0) Then
                  Write(LuPri,'(A)')
     &            'This seems to be an integer overflow!'
                  Call Cho_RWord2Byte(lX,Byte,Unt)
                  Write(LuPri,'(A,1P,D15.6,A,D15.6,1X,A,A)')
     &            'In double precision, lX=',lX,
     &            ' words (',Byte,Unt,')'
               End If
               irc=1
               lAvail=0
            Else
               lAvail=lWork-l_X
            End If
         Else
            lAvail = lWork - nT1amx ! all minus one vector (for reading)
         End If
         Call GAiGOp_Scal(lAvail,'min')
*        The argument ip_LnPQprod is only used for the case where full
*        Lpq-vectors are transformed for densities. Will be a dummy arg
*        for regular MP2.
         Accepted = ChoMP2_Setup_MemChk(iWork(ip_LnT1am),
     &                                  iWork(ip_LnPQprod),
     &                                  NumVec,nFrac,
     &                                  nSym,nBatch,lAvail)

         If (ForceBatch .and. nBatch.eq.1) Then
            Accepted = .false. ! force batching
         End If

         If (.not. Accepted) Then
            Call ChoMP2_deallocate(irc)
         End If

         If (irc.eq.1) Then
            Accepted=.false.
            nBatch=mBatch ! break while loop
         End If

      End Do

      If (.not. Accepted) Then
         Write(6,*) 'Accepted = false'
         irc = -1
         Return
      End If

C     Initialize file units.
C     ----------------------

      Do iSym = 1,nSym
         Do iTyp = 1,nTypF
            Call ChoMP2_OpenF(0,iTyp,iSym)
         End Do
      End Do

      End
*
      SubRoutine ChoMP2_Setup_Index(iFirst,iFirstS,NumOcc,
     &                              LnOcc,NumBatOrb,LnBatOrb,
     &                              LnT1am,LiT1am,
     &                              LnPQprod,LiPQprod,
     &                              LnMatij,LiMatij,
     &                              mSym,mBatch)
C
C     Thomas Bondo Pedersen, Nov. 2004 / Feb. 2005.
C
C     Purpose: set local index arrays and counters.
C
#include "implicit.fh"
      Integer iFirst(mBatch), NumOcc(mBatch)
      Integer iFirstS(mSym,mBatch), LnOcc(mSym,mBatch)
      Integer LnT1am(mSym,mBatch), LiT1am(mSym,mSym,mBatch)
      Integer LnMatij(mSym,mBatch), LiMatij(mSym,mSym,mBatch)
      Integer NumBatOrb(mBatch), LnBatOrb(mSym,mBatch)
      Integer LnPQprod(mSym,mBatch), LiPQprod(mSym,mSym,mBatch)
#include "cholesky.fh"
#include "chomp2_cfg.fh"
#include "chomp2.fh"

      Integer  Cho_iRange
      External Cho_iRange

      MulD2h(i,j)=iEor(i-1,j-1)+1

      If (mBatch .ne. nBatch) Then
         Call ChoMP2_Quit('ChoMP2_Setup_Index','mBatch !=  nBatch',
     &                    'Error')
      End If
      If (mSym .ne. nSym) Then
         Call ChoMP2_Quit('ChoMP2_Setup_Index','mSym !=  nSym',
     &                    'Error')
      End If

      Call Cho_iZero(iFirst,nBatch)
      Call Cho_iZero(iFirstS,nSym*nBatch)
      Call Cho_iZero(NumOcc,nBatch)
      Call Cho_iZero(NumBatOrb,nBatch)
      Call Cho_iZero(LnOcc,nSym*nBatch)
      Call Cho_iZero(LnBatOrb,nSym*nBatch)
      Call Cho_iZero(LnT1am,nSym*nBatch)
      Call Cho_iZero(LiT1am,nSym*nSym*nBatch)
      If(.false.) Then
         Call Cho_iZero(LnPQprod,nSym*nBatch)
         Call Cho_iZero(LiPQprod,nSym*nBatch)
      End If
      If (ChoAlg .eq. 2) Then
         Call Cho_iZero(LnMatij,nSym*nBatch)
         Call Cho_iZero(LiMatij,nSym*nSym*nBatch)
      End If
*
      Num = nBatOrbT/nBatch
*
*     I am not sure if NumOcc is used somewhere else so I will
*     define it as before even if Im using NumInBat for setting up
*     indeces in this routine. //Jonas
      Do iBatch = 1,nBatch
         If(.false.) Then
            NumBatOrb(iBatch) = Num
         Else
            NumOcc(iBatch) = Num
            NumBatOrb(iBatch) = Num
         End If
      End Do
*
      Left = nBatOrbT - nBatch*Num
      Do iBatch = nBatch,nBatch-Left+1,-1
         If(.false.) Then
            NumBatOrb(iBatch) = NumBatOrb(iBatch) + 1
         Else
            NumOcc(iBatch) = NumOcc(iBatch) + 1
            NumBatOrb(iBatch) = NumBatOrb(iBatch) + 1
         End If
      End Do
*
      iFirst(1) = 1
      Do i = 1,NumBatOrb(1)
         iSym = Cho_iRange(i,iBatOrb,nSym,.false.)
         If(.false.) Then
            LnBatOrb(iSym,1) = LnBatOrb(iSym,1) + 1
         Else
            LnBatOrb(iSym,1) = LnBatOrb(iSym,1) + 1
            LnOcc(iSym,1) = LnOcc(iSym,1) + 1
         End If
         If (iFirstS(iSym,1) .lt. 1) Then
            iFirstS(iSym,1) = i - iBatOrb(iSym)
         End If
      End Do

      Do iBatch = 2,nBatch
         iFirst(iBatch) = iFirst(iBatch-1) + NumBatOrb(iBatch-1)
         Do i = iFirst(iBatch),iFirst(iBatch)+NumBatOrb(iBatch)-1
            iSym = Cho_iRange(i,iBatOrb,nSym,.false.)
            If(.false.) Then
               LnBatOrb(iSym,iBatch) = LnBatOrb(iSym,iBatch) + 1
            Else
               LnBatOrb(iSym,iBatch) = LnBatOrb(iSym,iBatch) + 1
               LnOcc(iSym,iBatch) = LnOcc(iSym,iBatch) + 1
            End If
            If (iFirstS(iSym,iBatch) .lt. 1) Then
               iFirstS(iSym,iBatch) = i - iBatOrb(iSym)
            End If
         End Do
      End Do
*
      Do iBatch = 1,nBatch
         Do iSym = 1,nSym
            Do iSymi = 1,nSym
               iSyma = MulD2h(iSymi,iSym)
               LiT1am(iSyma,iSymi,iBatch) = LnT1am(iSym,iBatch)
               LnT1am(iSym,iBatch) = LnT1am(iSym,iBatch)
     &                             + nVir(iSyma)*LnOcc(iSymi,iBatch)
               If(.false.) Then
                  LiPQprod(iSyma,iSymi,iBatch) = LnPQprod(iSym,iBatch)
                  LnPQprod(iSym,iBatch) = LnPQprod(iSym,iBatch)
     &                                  + (nOcc(iSymA) + nVir(iSymA)
     &                                  +  nFro(iSymA) + nDel(iSymA))
     &                                  * (LnBatOrb(iSymI,iBatch))
               End If
            End Do
         End Do
      End Do

*
      If (ChoAlg .eq. 2) Then
         Do iBatch = 1,nBatch
            Do iSym = 1,nSym
               Do iSymj = 1,nSym
                  iSymi = MulD2h(iSymj,iSym)
                  If (iSymi .eq. iSymj) Then
                     LiMatij(iSymi,iSymi,iBatch) = LnMatij(iSym,iBatch)
                     LnMatij(iSym,iBatch) = LnMatij(iSym,iBatch)
     &                   + LnOcc(iSymi,iBatch)*(LnOcc(iSymi,iBatch)+1)/2
                  Else If (iSymi .lt. iSymj) Then
                     LiMatij(iSymi,iSymj,iBatch) = LnMatij(iSym,iBatch)
                     LiMatij(iSymj,iSymi,iBatch) = LnMatij(iSym,iBatch)
                     LnMatij(iSym,iBatch) = LnMatij(iSym,iBatch)
     &                   + LnOcc(iSymi,iBatch)*LnOcc(iSymj,iBatch)
                  End If
               End Do
            End Do
         End Do
      End If

      End
      Logical Function ChoMP2_Setup_MemChk(LnT1am,LnPQprod,NumVec,nFrac,
     &                                     nSym,nBatch,Mem)
C
C     Thomas Bondo Pedersen, Nov. 2004.
C
C     Purpose: Check memory availability.
C
#include "implicit.fh"
#include "chomp2_cfg.fh"
      Integer LnT1am(nSym,nBatch)
      Integer LnPQprod(nSym,nBatch)
      Integer NumVec(nSym), nFrac(nSym)

      Integer LnT2am, LiT2am(8)
      Integer LnPQRSprod,LiPQRSprod(8)

      Logical Accepted

      If (Mem .lt. 1) Then
         Accepted = .false.
         Go To 1 ! exit
      Else
         Accepted = .true.
      End If

      If (Laplace.and.SOS_MP2) Then
         xMem=dble(mem)
         xNeed=0.0d0
         Do iBatch=1,nBatch
            Do iSym=1,nSym
               Nai=LnT1am(iSym,iBatch)
               If (Nai.gt.0 .and. NumVec(iSym).gt.0) Then
                  xNeed=max(xNeed,dble(Nai)*dble(NumVec(iSym)))
               End If
            End Do
         End Do
         xLeft=xMem-xNeed
         If (xLeft.lt.1.0d0) Then
            Accepted=.false.
            Go To 1 ! exit
         End If
      Else
         Do iSym = 1,nSym
            If (nFrac(iSym) .lt. 1) Then
               Accepted = .false.
               Go To 1 ! exit
            End If
         End Do
         xMem = dble(mem)
         Do jBatch = 1,nBatch
            Do iBatch = 1,jBatch
               Call ChoMP2_Energy_GetInd(LnT2am,LiT2am,iBatch,jBatch)
               If(.false.) Then
                  Call ChoMP2_Energy_GetPQInd(LnPQRSprod,LiPQRSprod,
     &                                        iBatch,jBatch)
               End If
               If(.false.) Then
                  xInt  = dble(LnPQRSprod)
               Else
                  xInt  = dble(LnT2am)
               End If
               xLeft = xMem - xInt
               If (xInt.lt.1.0D0 .or. xLeft.lt.1.0D0) Then
                  Accepted = .false.
                  Go To 1 ! exit
               End If
               Do iSym = 1,nSym
                  If (iBatch .eq. jBatch) Then
                     If(.false.) Then
                        xDim = dble(LnPQprod(iSym,iBatch))
                     Else
                        xDim = dble(LnT1am(iSym,iBatch))
                     End If
                  Else
                     If(.false.) Then
                        xDim = dble(LnPQprod(iSym,iBatch))
     &                       + dble(LnPQprod(iSym,iBatch))
                     Else
                        xDim = dble(LnT1am(iSym,iBatch))
     &                       + dble(LnT1am(iSym,jBatch))
                     End If
                  End If
                  If (nFrac(iSym) .gt. NumVec(iSym)) Then
                     NumV = min(NumVec(iSym),1)
                  Else
                     NumV = NumVec(iSym)/nFrac(iSym)
                  End If
                  xNeed = xDim*dble(NumV)
                  xDiff = xLeft - xNeed
                  If (xDiff .lt. 1.0D0) Then
                     Accepted = .false.
                     Go To 1 ! exit
                  End If
               End Do
            End Do
         End Do
      End If

    1 ChoMP2_Setup_MemChk = Accepted
      End
      SubRoutine ChoMP2_Setup_Prt(irc)
C
C     Thomas Bondo Pedersen, Nov. 2004 / Feb. 2005.
C
C     Purpose: print setup for Cholesky MP2.
C
#include "implicit.fh"
#include "cholesky.fh"
#include "chomp2_cfg.fh"
#include "chomp2.fh"
#include "WrkSpc.fh"

      Integer iCount(8)

      iFirst(i)=iWork(ip_First-1+i)
      NumOcc(i)=iWork(ip_NumOcc-1+i)
      LnOcc(i,j)=iWork(ip_LnOcc-1+nSym*(j-1)+i)
      NumBatOrb(i)=iWork(ip_NumBatOrb-1+i)
      LnBatOrb(i,j)=iWork(ip_LnBatOrb-1+nSym*(j-1)+i)

      irc = 0

      Call Cho_iZero(iCount,nSym)

      Call Cho_Head('Cholesky MP2 Setup','=',80,6)
*     The values but not the names 'occupied' are updated to work
*     also for batching over all orbitals
      If (nBatch .gt. 1) Then
         Write(6,'(A,I6,A,I6,A)')
     &   'The list of',nBatOrbT,' occupied orbitals has been split in',
     &   nBatch,' batches:'
      Else If (nBatch .eq. 1) Then
         Write(6,'(A,I6,A)')
     &   'The list of',nBatOrbT,' occupied orbitals is not split:'
      Else
         Write(6,*) 'Oops, #batches over occupied orbitals ',
     &              'is non-positive: ',nBatch
         irc = -101
         Return
      End If

      Write(6,'(/,3X,A)')
     & ' Batch  First   Last #Occ/irrep'
      If (nSym .eq. 1) Then
         Write(6,'(3X,A)')
     &   '-------------------------------'
      Else If (nSym .eq. 2) Then
         Write(6,'(3X,A)')
     &   '-----------------------------------'
      Else If (nSym .eq. 4) Then
         Write(6,'(3X,A)')
     &   '-------------------------------------------------'
      Else If (nSym .eq. 8) Then
         Write(6,'(3X,A,A)')
     &   '------------------------------------------------------------',
     &   '-----------------'
      Else
         Write(6,*) 'Oops, #irreps is out of bounds: ',nSym
         irc = -102
         Return
      End If
      Do iBatch = 1,nBatch
         If(.false.) Then
            Write(6,'(3X,I6,1X,I6,1X,I6,8(1X,I6))')
     &           iBatch,iFirst(iBatch),
     &           iFirst(iBatch)+NumBatOrb(iBatch)-1,
     &           (LnBatOrb(iSym,iBatch),iSym=1,nSym)
         Else
            Write(6,'(3X,I6,1X,I6,1X,I6,8(1X,I6))')
     &           iBatch,iFirst(iBatch),iFirst(iBatch)+NumOcc(iBatch)-1,
     &           (LnOcc(iSym,iBatch),iSym=1,nSym)
         End If
         Do iSym = 1,nSym
            If(.false.) Then
               iCount(iSym) = iCount(iSym) + LnBatOrb(iSym,iBatch)
            Else
               iCount(iSym) = iCount(iSym) + LnOcc(iSym,iBatch)
            End If
         End Do
      End Do
      If (nSym .eq. 1) Then
         Write(6,'(3X,A)')
     &   '-------------------------------'
      Else If (nSym .eq. 2) Then
         Write(6,'(3X,A)')
     &   '-----------------------------------'
      Else If (nSym .eq. 4) Then
         Write(6,'(3X,A)')
     &   '-------------------------------------------------'
      Else If (nSym .eq. 8) Then
         Write(6,'(3X,A,A)')
     &   '------------------------------------------------------------',
     &   '-----------------'
      End If
      Write(6,'(3X,A,14X,8(1X,I6))') 'Total:',(iCount(iSym),iSym=1,nSym)
      If (nSym .eq. 1) Then
         Write(6,'(3X,A)')
     &   '-------------------------------'
      Else If (nSym .eq. 2) Then
         Write(6,'(3X,A)')
     &   '-----------------------------------'
      Else If (nSym .eq. 4) Then
         Write(6,'(3X,A)')
     &   '-------------------------------------------------'
      Else If (nSym .eq. 8) Then
         Write(6,'(3X,A,A)')
     &   '------------------------------------------------------------',
     &   '-----------------'
      End If
      Do iSym = 1,nSym
         If(.false.) Then
            If (iCount(iSym) .ne. nOcc(iSym)+nVir(iSym)
     &                          + nFro(iSym)+nDel(iSym)) Then
               Write(6,*) 'Oops, #Occ/irrep is incorrect....'
               irc = -103
               Return
            End If
         Else
            If (iCount(iSym) .ne. nOcc(iSym)) Then
               Write(6,*) 'Oops, #Occ/irrep is incorrect....'
               irc = -103
               Return
            End If
         End If
      End Do
      If (nBatch.gt.1 .and. ForceBatch) Then
         Write(6,'(/,A)')
     &   'Notice: batching has been requested by user.'
      End If

      Write(6,'(//,A)')
     & 'The following tasks will be performed:'
      Write(6,'(A)')
     & ' * AO-to-MO transformation of original Cholesky vectors.'
      If (DecoMP2) Then
         Write(6,'(A)')
     &   ' * Cholesky decomposition of (ai|bj) integrals.'
      End If
      If (nBatch .gt. 1) Then
         If (DecoMP2) Then
            Write(6,'(A)')
     &      ' * Presort of Cholesky vectors from (ai|bj) decomposition.'
         Else
            Write(6,'(A)') ' * Presort of MO Cholesky vectors.'
         End If
      End If
      If (Laplace.and.SOS_MP2) Then
         Write(6,'(A)')
     &   ' * Calculation of Laplace-SOS-MP2 correlation energy.'
         If (Laplace_nGridPoints.eq.0) Then
            Write(6,'(A)')
     &      '   Numerical Laplace integration quadrature: default'
         Else
            Write(6,'(A,I8)')
     &      '   Numerical Laplace integration quadrature:',
     &      Laplace_nGridPoints
         End If
      Else
         Write(6,'(A,A)')
     &   ' * On-the-fly assembly of (ai|bj) integrals and calculation ',
     &   'of MP2 energy correction.'
         Write(6,'(A,I3,A)')
     &   '   [Cholesky algorithm:',ChoAlg,']'
      End If

      Call xFlush(6)

      End
