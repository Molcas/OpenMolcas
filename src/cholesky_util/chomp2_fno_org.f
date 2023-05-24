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
! Copyright (C) 2008, Francesco Aquilante                              *
!***********************************************************************
      SubRoutine ChoMP2_fno_Org(irc,Delete,P_ab,P_ii,EOcc,EVir,Wrk,lWrk)
!
!      F. Aquilante, Geneva May 2008  (snick to Pedersen's code)
!
!
      use ChoMP2, only: iFirstS, LnOcc, LnT1am, LiMatij
      Implicit Real*8 (a-h,o-z)
      Logical Delete
      Real*8  P_ab(*), P_ii(*)
      Real*8  EOcc(*), EVir(*), Wrk(lWrk)
#include "cholesky.fh"
#include "chomp2_cfg.fh"
#include "chomp2.fh"
#include "chfnopt.fh"

      Character*7  ThisNm
      Character*14 SecNam
      Parameter (SecNam = 'ChoMP2_fno_Org', ThisNm = 'fno_Org')

      Integer nEnrVec(8), LnT2am, LiT2am(8), kP(8), lP(8)
      Integer nVaJi, iVaJi(8)

      Integer iDummy
      Parameter (iDummy = -999999)
      Real*8  xsDnom

      MulD2h(i,j)=iEor(i-1,j-1)+1
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j

      irc = 0

      kP(1)=1
      lP(1)=0
      Do iS=2,nSym
         kP(iS)=kP(iS-1)+nVir(iS-1)**2
         lP(iS)=lP(iS-1)+nOcc(iS-1)
      End Do

!     Set number and type of vectors.
!     -------------------------------

      If (DecoMP2) Then
         iTyp = 2
         Call iCopy(nSym,nMP2Vec,1,nEnrVec,1)
      Else
         iTyp = 1
         Call iCopy(nSym,NumCho,1,nEnrVec,1)
      End If

      If (MP2_small) Then

!        Loop over occupied orbital batches.
!        -----------------------------------

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
               Call FZero(Wrk(kMabij),LnT2am) ! initialize

!              Loop over Cholesky vector symmetries.
!              -------------------------------------

               Do iSym = 1,nSym

                  Nai = LnT1am(iSym,iBatch)
                  If (Nai.gt.0 .and. nEnrVec(iSym).gt.0) Then

!                    Reserve memory for reading a single vector.
!                    -------------------------------------------

                     kVecai = kEnd0
                     kEnd1  = kVecai + nT1am(iSym)
                     lWrk1  = lWrk   - kEnd1 + 1

                     If (lWrk1 .lt. Nai) Then
                        Call ChoMP2_Quit(SecNam,'Insufficient memory',
     &                                   '[ChoAlg.2.1]')
                     End If

!                    Set up batch over Cholesky vectors.
!                    -----------------------------------

                     nVec = min(lWrk1/Nai,nEnrVec(iSym))
                     If (nVec .lt. 1) Then ! should not happen
                        Call ChoMP2_Quit(SecNam,'Insufficient memory',
     &                                   '[ChoAlg.2.2]')
                     End If
                     nBat = (nEnrVec(iSym)-1)/nVec + 1

!                    Open Cholesky vector files.
!                    ---------------------------

                     Call ChoMP2_OpenF(1,iTyp,iSym)

!                    Start vector batch loop.
!                    ------------------------

                     Do iBat = 1,nBat

                        If (iBat .eq. nBat) Then
                           NumVec = nEnrVec(iSym) - nVec*(nBat-1)
                        Else
                           NumVec = nVec
                        End If
                        iVec1 = nVec*(iBat-1) + 1

!                       Set up index arrays for reordered vectors.
!                       ------------------------------------------

                        nVaJi = 0
                        Do iSymi = 1,nSym
                           iSyma = MulD2h(iSymi,iSym)
                           iVaJi(iSymi) = nVaJi
                           nVaJi = nVaJi
     &                          + nVir(iSyma)*NumVec*LnOcc(iSymi,iBatch)
                        End Do

!                       Pointer to reordered vectors: kVec.
!                       -----------------------------------

                        kVec  = kEnd1
                        kEnd2 = kVec  + nVaJi
                        lWrk2 = lWrk  - kEnd2 + 1
                        If (lWrk2 .lt. 0) Then ! should not happen
                           Call ChoMP2_Quit(SecNam,
     &                                      'Insufficient memory',
     &                                      '[ChoAlg.2.3]')
                        End If

!                       Read one vector at a time and reorder.
!                       --------------------------------------

                        iVec0 = iVec1 - 1
                        Do iVec = 1,NumVec

                           iOpt = 2
                           lTot = nT1am(iSym)
                           iAdr = nT1am(iSym)*(iVec0+iVec-1) + 1
                           Call ddaFile(lUnit_F(iSym,iTyp),iOpt,
     &                                  Wrk(kVecai),lTot,iAdr)

                           Do iSymi = 1,nSym
                              iSyma = MulD2h(iSymi,iSym)
                              i0 = iFirstS(iSymi,iBatch) - 1
                              Do i = 1,LnOcc(iSymi,iBatch)
                                 kOff1 = kVecai
     &                                 + iT1am(iSyma,iSymi)
     &                                 + nVir(iSyma)*(i0+i-1)
                                 kOff2 = kVec + iVaJi(iSymi)
     &                                 + nVir(iSyma)*NumVec*(i-1)
     &                                 + nVir(iSyma)*(iVec-1)
                                 Call dCopy_(nVir(iSyma),Wrk(kOff1),1,
     &                                      Wrk(kOff2),1)
                              End Do
                           End Do

                        End Do

!                       Compute M(ab,ii) .
!                       ---------------------------------------

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

!                    Close Cholesky vector files.
!                    ----------------------------

                     Call ChoMP2_OpenF(2,iTyp,iSym)


!                    Compute energy contribution
!                    ---------------------------------------

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

!                             Compute Energy contribution
!                             -------------------------------------
                              Do jb=1,nVir(iSymb)
                                 Do ja=1,nVir(iSymb)
                                    Dnom = EVir(iVir(iSymb)+ja)
     &                                   + EVir(iVir(iSymb)+jb)
     &                                   - 2.0d0*EOcc(iOcc(iSymj)+j)
                                    xsDnom = Dnom/(Dnom**2+shf**2)
                                    kOffMM = kOffM
     &                                     + nVir(iSymb)*(jb-1) +ja-1
                                    DeMP2 = DeMP2
!     &                                    + Wrk(kOffMM)**2/Dnom
     &                                    + Wrk(kOffMM)**2*xsDnom
                                 End Do
                              End Do

                           End Do

                        End If

                     End Do

                  End If

               End Do ! iSym

            End If

         End Do ! iBatch

!        If requested, delete vector files.
!        ----------------------------------

         If (Delete) Then
            Do iSym = 1,nSym
               Call ChoMP2_OpenF(1,iTyp,iSym)
               Call ChoMP2_OpenF(3,iTyp,iSym)
            End Do
         End If

         Return

      EndIf   !  MP2_small run

!     Loop over occupied orbital batches.
!     -----------------------------------

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
            Call FZero(Wrk(kMabij),LnT2am) ! initialize

!           Loop over Cholesky vector symmetries.
!           -------------------------------------

            Do iSym = 1,nSym

               Nai = LnT1am(iSym,iBatch)
               If (Nai.gt.0 .and. nEnrVec(iSym).gt.0) Then

!                 Reserve memory for reading a single vector.
!                 -------------------------------------------

                  kVecai = kEnd0
                  kEnd1  = kVecai + nT1am(iSym)
                  lWrk1  = lWrk   - kEnd1 + 1

                  If (lWrk1 .lt. Nai) Then
                     Call ChoMP2_Quit(SecNam,'Insufficient memory',
     &                                '[ChoAlg.2.1]')
                  End If

!                 Set up batch over Cholesky vectors.
!                 -----------------------------------

                  nVec = min(lWrk1/Nai,nEnrVec(iSym))
                  If (nVec .lt. 1) Then ! should not happen
                     Call ChoMP2_Quit(SecNam,'Insufficient memory',
     &                                '[ChoAlg.2.2]')
                  End If
                  nBat = (nEnrVec(iSym)-1)/nVec + 1

!                 Open Cholesky vector files.
!                 ---------------------------

                  Call ChoMP2_OpenF(1,iTyp,iSym)

!                 Start vector batch loop.
!                 ------------------------

                  Do iBat = 1,nBat

                     If (iBat .eq. nBat) Then
                        NumVec = nEnrVec(iSym) - nVec*(nBat-1)
                     Else
                        NumVec = nVec
                     End If
                     iVec1 = nVec*(iBat-1) + 1

!                    Set up index arrays for reordered vectors.
!                    ------------------------------------------

                     nVaJi = 0
                     Do iSymi = 1,nSym
                        iSyma = MulD2h(iSymi,iSym)
                        iVaJi(iSymi) = nVaJi
                        nVaJi = nVaJi
     &                       + nVir(iSyma)*NumVec*LnOcc(iSymi,iBatch)
                     End Do

!                    Pointer to reordered vectors: kVec.
!                    -----------------------------------

                     kVec  = kEnd1
                     kEnd2 = kVec  + nVaJi
                     lWrk2 = lWrk  - kEnd2 + 1
                     If (lWrk2 .lt. 0) Then ! should not happen
                        Call ChoMP2_Quit(SecNam,
     &                                   'Insufficient memory',
     &                                   '[ChoAlg.2.3]')
                     End If

!                    Read one vector at a time and reorder.
!                    --------------------------------------

                     iVec0 = iVec1 - 1
                     Do iVec = 1,NumVec

                        iOpt = 2
                        lTot = nT1am(iSym)
                        iAdr = nT1am(iSym)*(iVec0+iVec-1) + 1
                        Call ddaFile(lUnit_F(iSym,iTyp),iOpt,
     &                               Wrk(kVecai),lTot,iAdr)

                        Do iSymi = 1,nSym
                           iSyma = MulD2h(iSymi,iSym)
                           i0 = iFirstS(iSymi,iBatch) - 1
                           Do i = 1,LnOcc(iSymi,iBatch)
                              kOff1 = kVecai
     &                              + iT1am(iSyma,iSymi)
     &                              + nVir(iSyma)*(i0+i-1)
                              kOff2 = kVec + iVaJi(iSymi)
     &                              + nVir(iSyma)*NumVec*(i-1)
     &                              + nVir(iSyma)*(iVec-1)
                              Call dCopy_(nVir(iSyma),Wrk(kOff1),1,
     &                                   Wrk(kOff2),1)
                           End Do
                        End Do

                     End Do

!                    Compute M(ab,ii) .
!                    ---------------------------------------

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

                  End Do

                  Call Cho_GAdGOp(Wrk(kMabij),LnT2am,'+')

!                 Close Cholesky vector files.
!                 ----------------------------

                  Call ChoMP2_OpenF(2,iTyp,iSym)


                  Do iSymj = 1,nSym

                     iSymb = MulD2h(iSymj,iSym)

                     If (nVir(iSymb) .gt. 0) Then

                        Do j = 1,LnOcc(iSymj,iBatch)

                           i = j

                           ij = LiMatij(iSymj,iSymj,iBatch)
     &                        + iTri(i,j)

                           kOffM = kMabij + LiT2am(1)
     &                           + nMatab(1)*(ij-1)
     &                           + iMatab(iSymb,iSymb)

!                          Compute T(a,b)[i] and energy contrib.
!                          -------------------------------------
                           Do jb=1,nVir(iSymb)
                              Do ja=1,nVir(iSymb)
                                 Dnom = EVir(iVir(iSymb)+ja)
     &                                + EVir(iVir(iSymb)+jb)
     &                                - 2.0d0*EOcc(iOcc(iSymj)+j)
                                 xsDnom = Dnom/(Dnom**2+shf**2)
                                 kOffMM = kOffM
     &                                  + nVir(iSymb)*(jb-1) +ja-1
                                 DeMP2 = DeMP2
!     &                                 + Wrk(kOffMM)**2/Dnom
     &                                 + Wrk(kOffMM)**2/Dnom
!                                 Wrk(kOffMM) = Wrk(kOffMM)/Dnom
                                 Wrk(kOffMM) = Wrk(kOffMM)*xsDnom
                              End Do
                           End Do

                           P_ii(lP(iSymj)+i) = P_ii(lP(iSymj)+i)
     &                                       + ddot_(nVir(iSymb)**2,
     &                                                Wrk(kOffM),1,
     &                                                Wrk(kOffM),1)

!                          Compute P(a,b) += sum_c T(a,c)*T(c,b)
!                          -------------------------------------
                           Call dGeMM_('N','N',nVir(iSymb),
     &                                nVir(iSymb),nVir(iSymb),
     &                          1.0d0,Wrk(kOffM),nVir(iSymb),
     &                                Wrk(kOffM),nVir(iSymb),
     &                          1.0D0,P_ab(kP(iSymb)),nVir(iSymb))

                        End Do

                     End If

                  End Do

               End If

            End Do ! iSym

         End If

      End Do ! iBatch

!     If requested, delete vector files.
!     ----------------------------------

      If (Delete) Then
         Do iSym = 1,nSym
            Call ChoMP2_OpenF(1,iTyp,iSym)
            Call ChoMP2_OpenF(3,iTyp,iSym)
         End Do
      End If

      End
