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
      SubRoutine ChoMP2_Energy_Srt(irc,Delete,EMP2,EOcc,EVir,Wrk,lWrk)
C
C     Thomas Bondo Pedersen, Dec. 2004 / Feb. 2005.
C
C     Purpose: compute MP2 energy contribution using presorted MO
C              Cholesky vectors on disk.
C
#include "implicit.fh"
      Logical Delete
      Real*8  EOcc(*), EVir(*), Wrk(lWrk)
#include "cholesky.fh"
#include "chomp2_cfg.fh"
#include "chomp2.fh"
#include "WrkSpc.fh"

      Character*10 ThisNm
      Character*17 SecNam
      Parameter (SecNam = 'ChoMP2_Energy_Srt', ThisNm = 'Energy_Srt')

      Integer nEnrVec(8), LnT2am, LiT2am(8)
      Integer nVaJi, iVaJi(8)
      Integer iDummy
      Parameter (iDummy = -999999)

      Real*8 X(0:1)
      Data X /0.0D0,1.0D0/

      lUnit(i,j)=iWork(ip_lUnit-1+nSym*(j-1)+i)
      LnT1am(i,j)=iWork(ip_LnT1am-1+nSym*(j-1)+i)
      LiT1am(i,j,k)=iWork(ip_LiT1am-1+nSym*nSym*(k-1)+nSym*(j-1)+i)
      LiMatij(i,j,k)=iWork(ip_LiMatij-1+nSym*nSym*(k-1)+nSym*(j-1)+i)
      LnOcc(i,j)=iWork(ip_LnOcc-1+nSym*(j-1)+i)
      MulD2h(i,j)=iEor(i-1,j-1)+1
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j

      irc = 0

C     Set number of vectors.
C     ----------------------

      If (DecoMP2) Then
         Call iCopy(nSym,nMP2Vec,1,nEnrVec,1)
      Else
         Call iCopy(nSym,NumCho,1,nEnrVec,1)
      End If

C     Initialize MP2 energy correction.
C     ---------------------------------

      EMP2 = 0.0D0

C     Print header of status table.
C     -----------------------------

      If (Verbose) Then
         Call ChoMP2_Energy_Prt(SecNam,0,iDummy)
      End If

C     Loop over occupied orbital batches.
C     Special handling of diagonal batches for ChoAlg=2.
C     --------------------------------------------------

      Do iBatch = 1,nBatch
         If (Verbose) Then
            Call ChoMP2_Energy_Prt(SecNam,1,iBatch)
         End If
         Do jBatch = iBatch,nBatch

            Call ChoMP2_Energy_GetInd(LnT2am,LiT2am,iBatch,jBatch)

            kXaibj = 1
            kEnd0  = kXaibj + LnT2am
            lWrk0  = lWrk   - kEnd0 + 1
            If (lWrk0 .lt. 1) Then
               Call ChoMP2_Quit(SecNam,'insufficient memory','[0]')
            End If

C           Special code for iBatch=jBatch and ChoAlg=2:
C           compute M(ab,ij) = (ai|bj) with i<=j using level 3 BLAS.
C           For ChoAlg=1: use strictly lower triangular storage (=>
C           level 2 BLAS).
C           --------------------------------------------------------

            If (jBatch.eq.iBatch .and. ChoAlg.eq.2) Then

               kMabij = kXaibj  ! rename pointer
               Call Cho_dZero(Wrk(kMabij),LnT2am) ! initialize

C              Loop over Cholesky vector symmetries.
C              -------------------------------------

               Do iSym = 1,nSym

                  Nai = LnT1am(iSym,iBatch)
                  If (Nai.gt.0 .and. nEnrVec(iSym).gt.0) Then

C                    Reserve memory for reading a single vector.
C                    -------------------------------------------

                     kVecai = kEnd0
                     kEnd1  = kVecai + Nai
                     lWrk1  = lWrk   - kEnd1 + 1

                     If (lWrk1 .lt. Nai) Then
                        Call ChoMP2_Quit(SecNam,'Insufficient memory',
     &                                   '[ChoAlg.2.1]')
                     End If

C                    Set up batch over Cholesky vectors.
C                    -----------------------------------

                     nVec = min(lWrk1/Nai,nEnrVec(iSym))
                     If (nVec .lt. 1) Then ! should not happen
                        Call ChoMP2_Quit(SecNam,'Insufficient memory',
     &                                   '[ChoAlg.2.2]')
                     End If
                     nBat = (nEnrVec(iSym)-1)/nVec + 1

C                    Open Cholesky vector files.
C                    ---------------------------

                     Call ChoMP2_OpenB(1,iSym,iBatch)

C                    Start vector batch loop.
C                    ------------------------

                     Do iBat = 1,nBat

                        If (iBat .eq. nBat) Then
                           NumVec = nEnrVec(iSym) - nVec*(nBat-1)
                        Else
                           NumVec = nVec
                        End If
                        iVec1 = nVec*(iBat-1) + 1

C                       Set up index arrays for reordered vectors.
C                       ------------------------------------------

                        nVaJi = 0
                        Do iSymi = 1,nSym
                           iSyma = MulD2h(iSymi,iSym)
                           iVaJi(iSymi) = nVaJi
                           nVaJi = nVaJi
     &                          + nVir(iSyma)*NumVec*LnOcc(iSymi,iBatch)
                        End Do

C                       Pointer to reordered vectors: kVec.
C                       -----------------------------------

                        kVec  = kEnd1
                        kEnd2 = kVec  + nVaJi
                        lWrk2 = lWrk  - kEnd2 + 1
                        If (lWrk2 .lt. 0) Then ! should not happen
                           Call ChoMP2_Quit(SecNam,
     &                                      'Insufficient memory',
     &                                      '[ChoAlg.2.3]')
                        End If

C                       Read one vector at a time and reorder.
C                       --------------------------------------

                        iVec0 = iVec1 - 1
                        Do iVec = 1,NumVec

                           iOpt = 2
                           lTot = Nai
                           iAdr = Nai*(iVec0+iVec-1) + 1
                           Call ddaFile(lUnit(iSym,iBatch),iOpt,
     &                                  Wrk(kVecai),lTot,iAdr)

                           Do iSymi = 1,nSym
                              iSyma = MulD2h(iSymi,iSym)
                              Do i = 1,LnOcc(iSymi,iBatch)
                                 kOff1 = kVecai
     &                                 + LiT1am(iSyma,iSymi,iBatch)
     &                                 + nVir(iSyma)*(i-1)
                                 kOff2 = kVec + iVaJi(iSymi)
     &                                 + nVir(iSyma)*NumVec*(i-1)
     &                                 + nVir(iSyma)*(iVec-1)
                                 Call dCopy_(nVir(iSyma),Wrk(kOff1),1,
     &                                      Wrk(kOff2),1)
                              End Do
                           End Do

                        End Do

C                       Compute M(ab,ij) for i<=j.
C                       First do iSymi=iSymj, then iSymi<iSymj.
C                       ---------------------------------------

                        Do iSymj = 1,nSym

                           iSymb = MulD2h(iSymj,iSym)

                           If (nVir(iSymb) .gt. 0) Then

                              Do j = 1,LnOcc(iSymj,iBatch)
                                 Do i = 1,j

                                    ij = LiMatij(iSymj,iSymj,iBatch)
     &                                 + iTri(i,j)

                                    kOffi = kVec + iVaJi(iSymj)
     &                                    + nVir(iSymb)*NumVec*(i-1)
                                    kOffj = kVec + iVaJi(iSymj)
     &                                    + nVir(iSymb)*NumVec*(j-1)
                                    kOffM = kMabij + LiT2am(1)
     &                                    + nMatab(1)*(ij-1)
     &                                    + iMatab(iSymb,iSymb)

                                    Call DGEMM_('N','T',
     &                                   nVir(iSymb),nVir(iSymb),NumVec,
     &                                   1.0d0,Wrk(kOffi),nVir(iSymb),
     &                                         Wrk(kOffj),nVir(iSymb),
     &                                   1.0D0,Wrk(kOffM),nVir(iSymb))

                                 End Do
                              End Do

                              Do iSymi = 1,iSymj-1

                                 iSyma  = MulD2h(iSymi,iSym)
                                 iSymab = MulD2h(iSyma,iSymb)
                                 iSymij = iSymab

                                 If (LnOcc(iSymi,iBatch).gt.0 .and.
     &                               nVir(iSyma).gt.0) Then

                                    Do j = 1,LnOcc(iSymj,iBatch)
                                       Do i = 1,LnOcc(iSymi,iBatch)

                                          ij =
     &                                    LiMatij(iSymi,iSymj,iBatch)
     &                                    + LnOcc(iSymi,iBatch)*(j-1)
     &                                    + i

                                          kOffi = kVec + iVaJi(iSymi)
     &                                    + nVir(iSyma)*NumVec*(i-1)
                                          kOffj = kVec + iVaJi(iSymj)
     &                                    + nVir(iSymb)*NumVec*(j-1)
                                          kOffM = kMabij
     &                                    + LiT2am(iSymij)
     &                                    + nMatab(iSymab)*(ij-1)
     &                                    + iMatab(iSyma,iSymb)

                                          Call DGEMM_('N','T',
     &                                   nVir(iSyma),nVir(iSymb),NumVec,
     &                                    1.0d0,Wrk(kOffi),nVir(iSyma),
     &                                          Wrk(kOffj),nVir(iSymb),
     &                                    1.0D0,Wrk(kOffM),nVir(iSyma))

                                       End Do
                                    End Do

                                 End If

                              End Do

                           End If

                        End Do

                     End Do

C                    Close Cholesky vector files.
C                    ----------------------------

                     Call ChoMP2_OpenB(2,iSym,iBatch)

                  End If

               End Do ! iSym

            Else ! level 2 BLAS for diagonal batches.

C              Loop over Cholesky vector symmetries.
C              -------------------------------------

               Do iSym = 1,nSym

                  Nai = LnT1am(iSym,iBatch)
                  Nbj = LnT1am(iSym,jBatch)
                  If (Nai.gt.0 .and. Nbj.gt.0 .and. nEnrVec(iSym).gt.0)
     &            Then

C                    Setup Cholesky vector batching.
C                    -------------------------------

                     If (jBatch .eq. iBatch) Then
                       MinMem = Nai
                     Else
                       MinMem = Nai + Nbj
                     End If
                     NumVec = min(lWrk0/MinMem,nEnrVec(iSym))
                     If (NumVec .lt. 1) Then
                        Call ChoMP2_Quit(SecNam,'insufficient memory',
     &                                   '[1]')
                     End If
                     nBat = (nEnrVec(iSym) - 1)/NumVec + 1

C                    Open Cholesky vector files.
C                    ---------------------------

                     Call ChoMP2_OpenB(1,iSym,iBatch)
                     If (jBatch .ne. iBatch) Then
                        Call ChoMP2_OpenB(1,iSym,jBatch)
                     End If

C                    Cholesky vector batch loop.
C                    ---------------------------

                     Do iBat = 1,nBat

                        If (iBat .eq. nBat) Then
                           NumV = nEnrVec(iSym) - NumVec*(nBat-1)
                        Else
                           NumV = NumVec
                        End If
                        iVec1 = NumVec*(iBat-1) + 1

                        kVai = kEnd0
                        kVbj = kVai + Nai*NumV
                        If (jBatch .eq. iBatch) Then
                           kEnd1 = kVbj
                           kVbj  = kVai
                        Else
                           kEnd1 = kVbj + Nbj*NumV
                        End If
                        lWrk1 = lWrk - kEnd1 + 1
                        If (lWrk1 .lt. 0) Then ! this would be a bug...
                           Call ChoMP2_Quit(SecNam,
     &                                      'insufficient memory','[2]')
                        End If

C                       Read vectors.
C                       -------------

                        iOpt = 2
                        lTot = Nai*NumV
                        iAdr = Nai*(iVec1-1) + 1
                        Call ddaFile(lUnit(iSym,iBatch),iOpt,
     &                               Wrk(kVai),lTot,iAdr)
                        If (jBatch .ne. iBatch) Then
                           iOpt = 2
                           lTot = Nbj*NumV
                           iAdr = Nbj*(iVec1-1) + 1
                           Call ddaFile(lUnit(iSym,jBatch),iOpt,
     &                                  Wrk(kVbj),lTot,iAdr)
                        End If

C                       Compute integral contribution.
C                       ------------------------------

                        Fac   = X(min((iBat-1),1))
                        kXint = kXaibj + LiT2am(iSym)
                        If (iBatch .eq. jBatch) Then
                           Call dGeMM_Tri('N','T',Nai,Nai,NumV,
     &                                1.0D0,Wrk(kVai),Nai,Wrk(kVai),Nai,
     &                                Fac,Wrk(kXint),Nai)
                        Else
                           Call DGEMM_('N','T',Nai,Nbj,NumV,
     &                                1.0D0,Wrk(kVai),Nai,Wrk(kVbj),Nbj,
     &                                Fac,Wrk(kXint),Nai)
                        End If

                     End Do ! Cholesky vector batch

C                    Close Cholesky vector files.
C                    ----------------------------

                     Call ChoMP2_OpenB(2,iSym,iBatch)
                     If (jBatch .ne. iBatch) Then
                        Call ChoMP2_OpenB(2,iSym,jBatch)
                     End If

                  End If

               End Do ! iSym

            End If

C           Compute contribution to MP2 energy correction.
C           ----------------------------------------------

            Call ChoMP2_Energy_Contr(EMP2,EOcc,EVir,Wrk(kXaibj),
     &                               LnT2am,LiT2am,iBatch,jBatch)

         End Do ! jBatch
         If (Verbose) Then
            Call ChoMP2_Energy_Prt(SecNam,2,iBatch)
         End If
      End Do ! iBatch

C     Finish table.
C     -------------

      If (Verbose) Then
         Call ChoMP2_Energy_Prt(SecNam,3,iBatch)
      End If

C     Delete files if requested.
C     --------------------------

      If (Delete) Then
         Do iBatch = 1,nBatch
            Do iSym = 1,nSym
               Call ChoMP2_OpenB(1,iSym,iBatch)
               Call ChoMP2_OpenB(3,iSym,iBatch)
            End Do
         End Do
      End If

C     Change sign on energy.
C     ----------------------

      EMP2 = -EMP2

      End
