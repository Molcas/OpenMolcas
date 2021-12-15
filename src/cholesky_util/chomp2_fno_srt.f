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
* Copyright (C) 2008, Francesco Aquilante                              *
************************************************************************
      SubRoutine ChoMP2_fno_Srt(irc,Delete,P_ab,P_ii,EOcc,EVir,Wrk,lWrk)
C
C      F. Aquilante, Geneva May 2008  (snick to Pedersen's code)
C
C
#include "implicit.fh"
      Logical Delete
      Real*8  EOcc(*), EVir(*), Wrk(lWrk), P_ab(*), P_ii(*)
#include "cholesky.fh"
#include "chomp2_cfg.fh"
#include "chomp2.fh"
#include "WrkSpc.fh"
      Real*8  DeMP2
      Logical MP2_small
      Common / ChFNOPT/ DeMP2, MP2_small

      Character*10 ThisNm
      Character*17 SecNam
      Parameter (SecNam = 'ChoMP2_fno_Srt', ThisNm = 'fno_Srt')

      Integer nEnrVec(8), LnT2am, LiT2am(8), kP(8), lP(8)
      Integer nVaJi, iVaJi(8)
      Integer iDummy
      Parameter (iDummy = -999999)

      lUnit(i,j)=iWork(ip_lUnit-1+nSym*(j-1)+i)
      LnT1am(i,j)=iWork(ip_LnT1am-1+nSym*(j-1)+i)
      LiT1am(i,j,k)=iWork(ip_LiT1am-1+nSym*nSym*(k-1)+nSym*(j-1)+i)
      LiMatij(i,j,k)=iWork(ip_LiMatij-1+nSym*nSym*(k-1)+nSym*(j-1)+i)
      LnOcc(i,j)=iWork(ip_LnOcc-1+nSym*(j-1)+i)
      MulD2h(i,j)=iEor(i-1,j-1)+1
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j

      irc = 0

      kP(1)=1
      lP(1)=0
      Do iS=2,nSym
         kP(iS)=kP(iS-1)+nVir(iS-1)**2
         lP(iS)=lP(iS-1)+nOcc(iS-1)
      End Do

C     Set number of vectors.
C     ----------------------

      If (DecoMP2) Then
         Call iCopy(nSym,nMP2Vec,1,nEnrVec,1)
      Else
         Call iCopy(nSym,NumCho,1,nEnrVec,1)
      End If

      If (MP2_small) Then

C        Loop over occupied orbital batches.
C        --------------------------------------------------

         Do iBatch = 1,nBatch

            jBatch = iBatch

            Call ChoMP2_Energy_GetInd(LnT2am,LiT2am,iBatch,jBatch)

            kXaibj = 1
            kEnd0  = kXaibj + LnT2am
            lWrk0  = lWrk   - kEnd0 + 1
            If (lWrk0 .lt. 1) Then
               Call ChoMP2_Quit(SecNam,'insufficient memory','[0]')
            End If

            If (ChoAlg.eq.2) Then

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
     &                                '[ChoAlg.2.2]')
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
     &                                   Wrk(kOff2),1)
                              End Do
                           End Do

                        End Do

C                       Compute M(ab,ii) .
C                       ---------------------------------------

                        Do iSymj = 1,nSym

                           iSymb = MulD2h(iSymj,iSym)

                           If (nVir(iSymb) .gt. 0) Then

                              Do j = 1,LnOcc(iSymj,iBatch)

                                 i = j

                                 ij = LiMatij(iSymj,iSymj,iBatch)
     &                              + iTri(i,j)

                                 kOffi = kVec + iVaJi(iSymj)
     &                                 + nVir(iSymb)*NumVec*(i-1)
                                 kOffj = kVec + iVaJi(iSymj)
     &                                 + nVir(iSymb)*NumVec*(j-1)
                                 kOffM = kMabij + LiT2am(1)
     &                                 + nMatab(1)*(ij-1)
     &                                 + iMatab(iSymb,iSymb)

                                 Call dGeMM_('N','T',
     &                                nVir(iSymb),nVir(iSymb),NumVec,
     &                                1.0d0,Wrk(kOffi),nVir(iSymb),
     &                                      Wrk(kOffj),nVir(iSymb),
     &                                1.0D0,Wrk(kOffM),nVir(iSymb))

                              End Do

                           End If

                        End Do

                     End Do

                     Call Cho_GAdGOp(Wrk(kMabij),LnT2am,'+')

C                    Close Cholesky vector files.
C                    ----------------------------

                     Call ChoMP2_OpenB(2,iSym,iBatch)

                     Do iSymj = 1,nSym

                        iSymb = MulD2h(iSymj,iSym)

                        If (nVir(iSymb) .gt. 0) Then

                           Do j = 1,LnOcc(iSymj,iBatch)

                              i = j

                              ij = LiMatij(iSymj,iSymj,iBatch)
     &                           + iTri(i,j)

                              kOffM = kMabij + LiT2am(1)
     &                              + nMatab(1)*(ij-1)
     &                              + iMatab(iSymb,iSymb)

C                             Compute energy contribution
C                             -------------------------------------
                              Do jb=1,nVir(iSymb)
                                 Do ja=1,nVir(iSymb)
                                    Dnom = EVir(iVir(iSymb)+ja)
     &                                   + EVir(iVir(iSymb)+jb)
     &                                   - 2.0d0*EOcc(iOcc(iSymj)+j)
                                    kOffMM = kOffM
     &                                     + nVir(iSymb)*(jb-1) +ja-1
                                    DeMP2 = DeMP2
     &                                    + Wrk(kOffMM)**2/Dnom
                                 End Do
                              End Do

                           End Do

                        End If

                     End Do

                  End If

               End Do ! iSym

            End If

         End Do ! iBatch

C        Delete files if requested.
C        --------------------------

         If (Delete) Then
            Do iBatch = 1,nBatch
               Do iSym = 1,nSym
                  Call ChoMP2_OpenB(1,iSym,iBatch)
                  Call ChoMP2_OpenB(3,iSym,iBatch)
               End Do
            End Do
         End If

         Return

      EndIf


C     Loop over occupied orbital batches.
C     --------------------------------------------------

      Do iBatch = 1,nBatch

         jBatch = iBatch

         Call ChoMP2_Energy_GetInd(LnT2am,LiT2am,iBatch,jBatch)

         kXaibj = 1
         kEnd0  = kXaibj + LnT2am
         lWrk0  = lWrk   - kEnd0 + 1
         If (lWrk0 .lt. 1) Then
            Call ChoMP2_Quit(SecNam,'insufficient memory','[0]')
         End If

         If (ChoAlg.eq.2) Then

            kMabij = kXaibj  ! rename pointer
            Call Cho_dZero(Wrk(kMabij),LnT2am) ! initialize

C           Loop over Cholesky vector symmetries.
C           -------------------------------------

            Do iSym = 1,nSym

               Nai = LnT1am(iSym,iBatch)
               If (Nai.gt.0 .and. nEnrVec(iSym).gt.0) Then

C                 Reserve memory for reading a single vector.
C                 -------------------------------------------

                  kVecai = kEnd0
                  kEnd1  = kVecai + Nai
                  lWrk1  = lWrk   - kEnd1 + 1

                  If (lWrk1 .lt. Nai) Then
                     Call ChoMP2_Quit(SecNam,'Insufficient memory',
     &                                '[ChoAlg.2.1]')
                  End If

C                 Set up batch over Cholesky vectors.
C                 -----------------------------------

                  nVec = min(lWrk1/Nai,nEnrVec(iSym))
                  If (nVec .lt. 1) Then ! should not happen
                     Call ChoMP2_Quit(SecNam,'Insufficient memory',
     &                                '[ChoAlg.2.2]')
                  End If
                  nBat = (nEnrVec(iSym)-1)/nVec + 1

C                 Open Cholesky vector files.
C                 ---------------------------

                  Call ChoMP2_OpenB(1,iSym,iBatch)

C                 Start vector batch loop.
C                 ------------------------

                  Do iBat = 1,nBat

                     If (iBat .eq. nBat) Then
                        NumVec = nEnrVec(iSym) - nVec*(nBat-1)
                     Else
                        NumVec = nVec
                     End If
                     iVec1 = nVec*(iBat-1) + 1

C                    Set up index arrays for reordered vectors.
C                    ------------------------------------------

                     nVaJi = 0
                     Do iSymi = 1,nSym
                        iSyma = MulD2h(iSymi,iSym)
                        iVaJi(iSymi) = nVaJi
                        nVaJi = nVaJi
     &                       + nVir(iSyma)*NumVec*LnOcc(iSymi,iBatch)
                     End Do

C                    Pointer to reordered vectors: kVec.
C                    -----------------------------------

                     kVec  = kEnd1
                     kEnd2 = kVec  + nVaJi
                     lWrk2 = lWrk  - kEnd2 + 1
                     If (lWrk2 .lt. 0) Then ! should not happen
                        Call ChoMP2_Quit(SecNam,
     &                                   'Insufficient memory',
     &                                   '[ChoAlg.2.3]')
                     End If

C                    Read one vector at a time and reorder.
C                    --------------------------------------

                     iVec0 = iVec1 - 1
                     Do iVec = 1,NumVec

                        iOpt = 2
                        lTot = Nai
                        iAdr = Nai*(iVec0+iVec-1) + 1
                        Call ddaFile(lUnit(iSym,iBatch),iOpt,
     &                               Wrk(kVecai),lTot,iAdr)

                        Do iSymi = 1,nSym
                           iSyma = MulD2h(iSymi,iSym)
                           Do i = 1,LnOcc(iSymi,iBatch)
                              kOff1 = kVecai
     &                              + LiT1am(iSyma,iSymi,iBatch)
     &                              + nVir(iSyma)*(i-1)
                              kOff2 = kVec + iVaJi(iSymi)
     &                              + nVir(iSyma)*NumVec*(i-1)
     &                              + nVir(iSyma)*(iVec-1)
                              Call dCopy_(nVir(iSyma),Wrk(kOff1),1,
     &                                   Wrk(kOff2),1)
                           End Do
                        End Do

                     End Do

C                    Compute M(ab,ii) .
C                    ---------------------------------------

                     Do iSymj = 1,nSym

                        iSymb = MulD2h(iSymj,iSym)

                        If (nVir(iSymb) .gt. 0) Then

                           Do j = 1,LnOcc(iSymj,iBatch)

                              i = j

                              ij = LiMatij(iSymj,iSymj,iBatch)
     &                           + iTri(i,j)

                              kOffi = kVec + iVaJi(iSymj)
     &                              + nVir(iSymb)*NumVec*(i-1)
                              kOffj = kVec + iVaJi(iSymj)
     &                              + nVir(iSymb)*NumVec*(j-1)
                              kOffM = kMabij + LiT2am(1)
     &                              + nMatab(1)*(ij-1)
     &                              + iMatab(iSymb,iSymb)

                              Call dGeMM_('N','T',
     &                             nVir(iSymb),nVir(iSymb),NumVec,
     &                             1.0d0,Wrk(kOffi),nVir(iSymb),
     &                                   Wrk(kOffj),nVir(iSymb),
     &                             1.0D0,Wrk(kOffM),nVir(iSymb))

                           End Do

                        End If

                     End Do

                     Call Cho_GAdGOp(Wrk(kMabij),LnT2am,'+')

C                    Close Cholesky vector files.
C                    ----------------------------

                     Call ChoMP2_OpenB(2,iSym,iBatch)

                     Do iSymj = 1,nSym

                        iSymb = MulD2h(iSymj,iSym)

                        If (nVir(iSymb) .gt. 0) Then

                           Do j = 1,LnOcc(iSymj,iBatch)

                              i = j

                              ij = LiMatij(iSymj,iSymj,iBatch)
     &                           + iTri(i,j)

                              kOffM = kMabij + LiT2am(1)
     &                              + nMatab(1)*(ij-1)
     &                              + iMatab(iSymb,iSymb)

C                             Compute T(a,b)[i] and energy contrib.
C                             -------------------------------------
                              Do jb=1,nVir(iSymb)
                                 Do ja=1,nVir(iSymb)
                                    Dnom = EVir(iVir(iSymb)+ja)
     &                                   + EVir(iVir(iSymb)+jb)
     &                                   - 2.0d0*EOcc(iOcc(iSymj)+j)
                                    kOffMM = kOffM
     &                                     + nVir(iSymb)*(jb-1) +ja-1
                                    DeMP2 = DeMP2
     &                                    + Wrk(kOffMM)**2/Dnom
                                    Wrk(kOffMM) = Wrk(kOffMM)/Dnom
                                 End Do
                              End Do

                              P_ii(lP(iSymj)+i) = P_ii(lP(iSymj)+i)
     &                                          + ddot_(nVir(iSymb)**2,
     &                                                   Wrk(kOffM),1,
     &                                                   Wrk(kOffM),1)

C                             Compute P(a,b) += sum_c T(a,c)*T(c,b)
C                             -------------------------------------
                              Call dGeMM_('N','N',nVir(iSymb),
     &                                   nVir(iSymb),nVir(iSymb),
     &                             1.0d0,Wrk(kOffM),nVir(iSymb),
     &                                   Wrk(kOffM),nVir(iSymb),
     &                             1.0D0,P_ab(kP(iSymb)),nVir(iSymb))

                           End Do

                        End If

                     End Do

                  End Do

               End If

            End Do ! iSym

         End If

      End Do ! iBatch

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

      End
