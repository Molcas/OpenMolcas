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
* Copyright (C) 2008, Jonas Bostrom                                    *
************************************************************************

      SubRoutine ChoMP2_Read_Batch(LnPQRSprod,LiPQRSprod,Wrk,lWrk,
     &                             iBatch,jBatch,kXpqrs)
C
C     Jonas Bostrom, Aug 2008. (Generalization and modification of some
C                               stuff in ChoMP2_energy_*)
C
C     Purpose: Reads parts of cholesky vectors from disk and
C              multiply them into two integral batches (or one batch if
C              jBatch=iBatch).
C
      Implicit Real*8 (a-h,o-z)
      Real*8 Wrk(lWrk)
#include "cholesky.fh"
#include "chomp2_cfg.fh"
#include "chomp2.fh"
#include "WrkSpc.fh"
*
      Character*15 ThisNm
      Character*25 SecNam
      Parameter (SecNam = 'ChoMP2_Read_Batch',
     &           ThisNm = 'Read_Batch')
*
      Integer nEnrVec(8), LnPQRSprod, LiPQRSprod(8)
      Real*8 X(0:1)
      Data X /0.0D0,1.0D0/
*
      LnPQprod(i,j)=iWork(ip_LnPQprod-1+nSym*(j-1)+i)
*
C     Set number and type of vectors.
C     -------------------------------
*
      iTyp = 1
      Call iCopy(nSym,NumCho,1,nEnrVec,1)
*
*     Allocate memory for integrals.
*     ------------------------------
*
      kXpqrs = 1
      kEnd0  = kXpqrs + LnPQRSprod
      lWrk0  = lWrk   - kEnd0 + 1
      If (lWrk0 .lt. 1) Then
         Call ChoMP2_Quit(SecNam,'insufficient memory','[0]')
      End If
*
C     Special code for iBatch=jBatch and ChoAlg=2:
C     compute M(ab,ij) = (ai|bj) with i<=j using level 3 BLAS.
C     For ChoAlg=1: use strictly lower triangular storage (=>
C     level 2 BLAS).
C     --------------------------------------------------------
*
      If (ChoAlg.eq.2) Then
         Write(6,*) 'No support for Cholesky algorithm 2'
      Else
*
C     Loop over Cholesky vector symmetries.
C     -------------------------------------
*
         Do iSym = 1,nSym
*
            Nai = LnPQprod(iSym,iBatch)
            Nbj = LnPQprod(iSym,jBatch)
            If ((Nai.gt.0) .and. (Nbj.gt.0) .and.
     &           (nEnrVec(iSym).gt.0)) Then
*
*
C     Allocate memory for reading 1 vector.
C     -------------------------------------
*
               kRead = kEnd0
               If(nBatch .ne. 1) Then
                  kEnd1 = kRead + nPQ_prod(iSym)
                  lWrk1 = lWrk  - kEnd1 + 1
                  If (lWrk1 .lt. 1) Then
                     Call ChoMP2_Quit(SecNam,'insufficient memory',
     &                    '[0.1]')
                  End If
               Else
                  kEnd1 = kRead
                  lWrk1 = lWrk0
               End If
*
C     Setup Cholesky vector batching.
C     -------------------------------
*
               If (jBatch .eq. iBatch) Then
                  MinMem = Nai
               Else
                  MinMem = Nai + Nbj
               End If
               NumVec = min(lWrk1/MinMem,nEnrVec(iSym))
               If (NumVec .lt. 1) Then
                  Call ChoMP2_Quit(SecNam,'insufficient memory',
     &                 '[1]')
               End If
               nBat = (nEnrVec(iSym) - 1)/NumVec + 1
*
*     Open Cholesky vector file.
*     --------------------------
*
               Call ChoMP2_OpenF(1,iTyp,iSym)
*
C     Cholesky vector batch loop.
C     ---------------------------
*
               Do iBat = 1,nBat
*
                  If (iBat .eq. nBat) Then
                     NumV = nEnrVec(iSym) - NumVec*(nBat-1)
                  Else
                     NumV = NumVec
                  End If
                  iVec1 = NumVec*(iBat-1) + 1
*
                  kVai = kEnd1
                  kVbj = 0
                  kEnd2 = 0
                  lWrk2 = 0
                  If(nBatch .ne. 1) Then
                     kVbj = kVai + Nai*NumV
                     If (jBatch .eq. iBatch) Then
                        kEnd2 = kVbj
                        kVbj  = kVai
                     Else
                        kEnd2 = kVbj + Nbj*NumV
                     End If
                     lWrk2 = lWrk - kEnd2 + 1
                     If (lWrk2 .lt. 0) Then ! this would be a bug...
                        Call ChoMP2_Quit(SecNam,
     &                       'insufficient memory',
     &                       '[2]')
                     End If
                  End If

*
C     Read vectors, copy out sub-blocks.
C     ----------------------------------
*
                  If(nBatch.eq.1) Then
*
                     iOpt = 2
                     lTot = nPQ_prod(iSym)*NumV
                     iAdr = nPQ_prod(iSym)*(iVec1 - 1) + 1
                     Call ddaFile(lUnit_F(iSym,iTyp),iOpt,Wrk(kVai),
     &                            lTot,iAdr)
                  Else
                     Do iVec = 1,NumV
*
                        jVec = iVec1 + iVec - 1
                        iOpt = 2
                        lTot = nPQ_prod(iSym)
                        iAdr = nPQ_prod(iSym)*(jVec-1) + 1
                        Call ddaFile(lUnit_F(iSym,iTyp),iOpt,
     &                       Wrk(kRead),lTot,iAdr)
                        kOff = kVai + Nai*(iVec - 1)
                        Call ChoMP2_Srt(Wrk(kRead),Wrk(kOff),
     &                                  1,iSym,iBatch)
                        If (jBatch .ne. iBatch) Then
                           kOff = kVbj + Nbj*(iVec - 1)
                           Call ChoMP2_Srt(Wrk(kRead),Wrk(kOff),
     &                          1,iSym,jBatch)
                        End If
                     End Do
                  End If
*
C     Compute integral contribution.
C     ------------------------------
*
                  Fac   = X(min((iBat-1),1))
                  kXint = kXpqrs + LiPQRSprod(iSym)
                  If (iBatch .eq. jBatch) Then
                     Call dGeMM_Tri('N','T',Nai,Nai,NumV,
     &                    1.0D0,Wrk(kVai),Nai,
     &                    Wrk(kVai),Nai,
     &                    Fac,Wrk(kXint),Nai)
                  Else
                     Call DGEMM_('N','T',Nai,Nbj,NumV,
     &                    1.0D0,Wrk(kVai),Nai,Wrk(kVbj),Nbj,
     &                    Fac,Wrk(kXint),Nai)
                  End If
               End Do           ! Cholesky vector batch
*
C     Close Cholesky vector files.
C     ----------------------------
*
               Call ChoMP2_OpenF(2,iTyp,iSym)
*
            End If
*
         End Do                 ! iSym
*
      End If
*
      Return
      End
