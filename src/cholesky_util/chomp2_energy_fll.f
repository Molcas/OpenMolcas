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
* Copyright (C) 2004, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine ChoMP2_Energy_Fll(irc,Delete,EMP2,EOcc,EVir,Wrk,lWrk)
C
C     Thomas Bondo Pedersen, Dec. 2004.
C
C     Purpose: compute MP2 energy contribution using original
C              Cholesky vectors on disk for the case nBatch=1.
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
      Parameter (SecNam = 'ChoMP2_Energy_Fll', ThisNm = 'Energy_Fll')

      Integer nEnrVec(8), LnT2am, LiT2am(8)
      Integer nVaJi, iVaJi(8)

      Real*8 X(0:1)
      Data X /0.0D0,1.0D0/

      LiMatij(i,j,k)=iWork(ip_LiMatij-1+nSym*nSym*(k-1)+nSym*(j-1)+i)
      MulD2h(i,j)=iEor(i-1,j-1)+1
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j

      If (nBatch .ne. 1) Then
         irc = -1
         Return
      End If

      Call qEnter(ThisNm)
      irc = 0

C     Determine if vector files are to be deleted after use.
C     ------------------------------------------------------

      If (Delete) Then
         iClos = 3
      Else
         iClos = 2
      End If

C     Set number and type of vectors.
C     -------------------------------

      If (DecoMP2) Then
         iTyp = 2
         Call iCopy(nSym,nMP2Vec,1,nEnrVec,1)
      Else
         iTyp = 1
         Call iCopy(nSym,NumCho,1,nEnrVec,1)
      End If

C     Set (ai|bj) indices.
C     --------------------

      Call ChoMP2_Energy_GetInd(LnT2am,LiT2am,1,1)

C     Allocate memory for integrals.
C     ------------------------------

      kXaibj = 1
      kEnd0  = kXaibj + LnT2am
      lWrk0  = lWrk   - kEnd0 + 1
      If (lWrk0 .lt. 0) Then
         Call ChoMP2_Quit(SecNam,'insufficient memory','[0]')
      End If

C     Initialize MP2 energy correction.
C     ---------------------------------

      EMP2 = 0.0D0

C     Special code for ChoAlg=2:
C     compute M(ab,ij) = (ai|bj) with i<=j using level 3 BLAS.
C     For ChoAlg=1: use strictly lower triangular storage (=>
C     level 2 BLAS).
C     --------------------------------------------------------

      If (ChoAlg .eq. 2) Then ! level 3 BLAS algorithm

         kMabij = kXaibj ! rename pointer
         Call Cho_dZero(Wrk(kMabij),LnT2am) ! initialize

C        Loop over Cholesky symmetries.
C        ------------------------------

         Do iSym = 1,nSym

            Nai = nT1am(iSym)
            If (Nai.gt.0 .and. nEnrVec(iSym).gt.0) Then

C              Reserve memory for reading a single vector.
C              -------------------------------------------

               kVecai = kEnd0
               kEnd1  = kVecai + Nai
               lWrk1  = lWrk   - kEnd1 + 1

               If (lWrk1 .lt. Nai) Then
                  Call ChoMP2_Quit(SecNam,'Insufficient memory',
     &                             '[ChoAlg.2.1]')
               End If

C              Set up batch over Cholesky vectors.
C              -----------------------------------

               nVec = min(lWrk1/Nai,nEnrVec(iSym))
               If (nVec .lt. 1) Then ! should not happen
                  Call ChoMP2_Quit(SecNam,'Insufficient memory',
     &                             '[ChoAlg.2.2]')
               End If
               nBat = (nEnrVec(iSym)-1)/nVec + 1

C              Open Cholesky vector files.
C              ---------------------------

               Call ChoMP2_OpenF(1,iTyp,iSym)

C              Start vector batch loop.
C              ------------------------

               Do iBat = 1,nBat

                  If (iBat .eq. nBat) Then
                     NumVec = nEnrVec(iSym) - nVec*(nBat-1)
                  Else
                     NumVec = nVec
                  End If
                  iVec1 = nVec*(iBat-1) + 1

C                 Set up index arrays for reordered vectors.
C                 ------------------------------------------

                  nVaJi = 0
                  Do iSymi = 1,nSym
                     iSyma = MulD2h(iSymi,iSym)
                     iVaJi(iSymi) = nVaJi
                     nVaJi = nVaJi
     &                    + nVir(iSyma)*NumVec*nOcc(iSymi)
                  End Do

C                 Pointer to reordered vectors: kVec.
C                 -----------------------------------

                  kVec  = kEnd1
                  kEnd2 = kVec  + nVaJi
                  lWrk2 = lWrk  - kEnd2 + 1
                  If (lWrk2 .lt. 0) Then ! should not happen
                     Call ChoMP2_Quit(SecNam,
     &                                'Insufficient memory',
     &                                '[ChoAlg.2.3]')
                  End If

C                 Read one vector at a time and reorder.
C                 --------------------------------------

                  iVec0 = iVec1 - 1
                  Do iVec = 1,NumVec

                     iOpt = 2
                     lTot = nT1am(iSym)
                     iAdr = nT1am(iSym)*(iVec0+iVec-1) + 1
                     Call ddaFile(lUnit_F(iSym,iTyp),iOpt,
     &                            Wrk(kVecai),lTot,iAdr)

                     Do iSymi = 1,nSym
                        iSyma = MulD2h(iSymi,iSym)
                        Do i = 1,nOcc(iSymi)
                           kOff1 = kVecai
     &                           + iT1am(iSyma,iSymi)
     &                           + nVir(iSyma)*(i-1)
                           kOff2 = kVec + iVaJi(iSymi)
     &                           + nVir(iSyma)*NumVec*(i-1)
     &                           + nVir(iSyma)*(iVec-1)
                           Call dCopy_(nVir(iSyma),Wrk(kOff1),1,
     &                                Wrk(kOff2),1)
                        End Do
                     End Do

                  End Do

C                 Compute M(ab,ij) for i<=j.
C                 First do iSymi=iSymj, then iSymi<iSymj.
C                 ---------------------------------------

                  Do iSymj = 1,nSym

                     iSymb = MulD2h(iSymj,iSym)

                     If (nVir(iSymb) .gt. 0) Then

                        Do j = 1,nOcc(iSymj)
                           Do i = 1,j

                              ij = LiMatij(iSymj,iSymj,1)
     &                           + iTri(i,j)

                              kOffi = kVec + iVaJi(iSymj)
     &                              + nVir(iSymb)*NumVec*(i-1)
                              kOffj = kVec + iVaJi(iSymj)
     &                              + nVir(iSymb)*NumVec*(j-1)
                              kOffM = kMabij + LiT2am(1)
     &                              + nMatab(1)*(ij-1)
     &                              + iMatab(iSymb,iSymb)

                              Call DGEMM_('N','T',
     &                             nVir(iSymb),nVir(iSymb),NumVec,
     &                             1.0d0,Wrk(kOffi),nVir(iSymb),
     &                                   Wrk(kOffj),nVir(iSymb),
     &                             1.0D0,Wrk(kOffM),nVir(iSymb))

                           End Do
                        End Do

                        Do iSymi = 1,iSymj-1

                           iSyma  = MulD2h(iSymi,iSym)
                           iSymab = MulD2h(iSyma,iSymb)
                           iSymij = iSymab

                           If (nOcc(iSymi).gt.0 .and.
     &                         nVir(iSyma).gt.0) Then

                              Do j = 1,nOcc(iSymj)
                                 Do i = 1,nOcc(iSymi)

                                    ij = LiMatij(iSymi,iSymj,1)
     &                                 + nOcc(iSymi)*(j-1) + i

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

C              Close (and possibly delete) Cholesky vector files.
C              --------------------------------------------------

               Call ChoMP2_OpenF(iClos,iTyp,iSym)

            End If

         End Do ! iSym

      Else ! level 2 BLAS algorithm

C        Loop over Cholesky symmetries.
C        ------------------------------

         Do iSym = 1,nSym
            If (nT1am(iSym).gt.0 .and. nEnrVec(iSym).gt.0) Then

C              Set up Cholesky batch.
C              ----------------------

               NumVec = min(lWrk0/nT1am(iSym),nEnrVec(iSym))
               If (NumVec .lt. 1) Then
                  Call ChoMP2_Quit(SecNam,'insufficient memory','[2]')
               End IF
               nBat = (nEnrVec(iSym) - 1)/NumVec + 1

C              Open Cholesky vector file.
C              --------------------------

               Call ChoMP2_OpenF(1,iTyp,iSym)

C              Start batch loop.
C              -----------------

               kOffX = kXaibj + LiT2am(iSym)
               Do iBat = 1,nBat

                  If (iBat .eq. nBat) Then
                     NumV = nEnrVec(iSym) - NumVec*(nBat-1)
                  Else
                     NumV = NumVec
                  End If
                  iVec1 = NumVec*(iBat-1) + 1

C                 Read vectors.
C                 -------------

                  kVec = kEnd0

                  iOpt = 2
                  lTot = nT1am(iSym)*NumV
                  iAdr = nT1am(iSym)*(iVec1 - 1) + 1
                  Call ddaFile(lUnit_F(iSym,iTyp),iOpt,Wrk(kVec),lTot,
     &                         iAdr)

C                 Compute contribution.
C                 ---------------------

                  Fac   = X(min((iBat-1),1))
                  kXint = kXaibj + LiT2am(iSym)
                  Call dGeMM_Tri('N','T',nT1am(iSym),nT1am(iSym),NumV,
     &                           1.0D0,Wrk(kVec),nT1am(iSym),
     &                           Wrk(kVec),nT1am(iSym),
     &                           Fac,Wrk(kXint),nT1am(iSym))

               End Do

C              Close (and possibly delete) Cholesky vector file.
C              -------------------------------------------------

               Call ChoMP2_OpenF(iClos,iTyp,iSym)

            End If
         End Do

      End If ! ChoAlg

C     Compute energy contribution.
C     ----------------------------

      Call ChoMP2_Energy_Contr(EMP2,EOcc,EVir,Wrk(kXaibj),
     &                         LnT2am,LiT2am,1,1)

C     Change sign on energy.
C     ----------------------

      EMP2 = -EMP2

      Call qExit(ThisNm)
      End
