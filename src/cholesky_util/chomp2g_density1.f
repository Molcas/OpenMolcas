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
* Copyright (C) 2010, Jonas Bostrom                                    *
************************************************************************
      SubRoutine ChoMP2g_density1(irc,EOcc,EVir,EFro,Wrk,lWrk)

*     Jonas Bostrom, Feb 2010
*
*     Purpose: To Compute MP2 density from Cholesky MO-vectors and
*              decomposed MP2 amplitudes.

#include "implicit.fh"
#include "chomp2g.fh"
#include "chomp2.fh"
#include "chomp2_cfg.fh"
#include "cholesky.fh"
#include "choorb.fh"
#include "WrkSpc.fh"
*
      Character*16 SecNam
      Parameter (SecNam='ChoMP2g_density1')
*
      Character Fname*5
      Real*8 Wrk(lWrk), EOcc(*), EVir(*), EFro(*)

      Real*8 X(0:1)
      Data X /0.0D0,1.0D0/
*     -----------------------------
      MulD2h(i,j)=iEor(i-1,j-1) + 1
      iAdrVec(i,j,k) = (i-1) + (j-1)*nSym + (k-1)*nSym*nSym
*     -----------------------------
      maxvalue = 200

*     Do not delete vectors
*     ---------------------
      iClos = 2

*     Set type of Choleskyvectors
*     ---------------------------
      iTypL = 1
      iTypR = 2
      iVecOV = 6
      iVecOF = 4
      iVecOO = 5
      iVecFV = 3
      iVecVV = 9
*
*     Calc max number of cholesky vectors in a specific sym
*     -----------------------------------------------------
      nMP2VecMax = 0
      NumChoMax = 0
      Do i = 1, nSym
         nMP2VecMax = Max(nMP2VecMax,nMP2Vec(i))
         NumChoMax = Max(NumChoMax, NumCho(i))
      End Do
*
      nMoMoMax = 0
      Do iSym = 1, nSym
         nMoMoMax = max(nMoMo(iSym,iVecOV),nMoMoMax)
      End Do

*     Allocate Memory for Pab-Vector
*     ------------------------------
      lPab = nVir(1)*nVir(1)
      kPab(1) = 1
      Do iSym = 2, nSym
         kPab(iSym) = kPab(iSym-1) + nVir(iSym-1)*nVir(iSym-1)
         lPab = lPab + nVir(iSym)*nVir(iSym)
      End Do
      kEndPab = kPab(1) + lPab
      Call FZero(Wrk(kPab(1)),lPab)

*     Allocate Memory for Wab-Vector
*     ------------------------------
      lWab = nVir(1)*nVir(1)
      kWab(1) = kEndPab
      Do iSym = 2, nSym
         kWab(iSym) = kWab(iSym-1) + nVir(iSym-1)*nVir(iSym-1)
         lWab = lWab + nVir(iSym)*nVir(iSym)
      End Do
      kEndWab = kWab(1) + lWab
      Call FZero(Wrk(kWab(1)),lWab)

*     Allocate memory for Pij-Vector
*     ------------------------------

      lPij = nOcc(1)*nOcc(1)
      kPij(1) = kEndWab
      Do iSym = 2, nSym
         kPij(iSym) = kPij(iSym-1) + nOcc(iSym-1)*nOcc(iSym-1)
         lPij = lPij + nOcc(iSym)*nOcc(iSym)
      End Do
      kEndPij = kPij(1) + lPij
      Call FZero(Wrk(kPij(1)),lPij)

*     Allocate memory for Wij-Vector
*     ------------------------------

      lWij = nOcc(1)*nOcc(1)
      kWij(1) = kEndPij
      Do iSym = 2, nSym
         kWij(iSym) = kWij(iSym-1) + nOcc(iSym-1)*nOcc(iSym-1)
         lWij = lWij + nOcc(iSym)*nOcc(iSym)
      End Do
      kEndWij = kWij(1) + lWij
      Call FZero(Wrk(kWij(1)),lWij)

*     Allocate memory for Pai-Vector
*     ------------------------------

      lPai = nVir(1)*nOcc(1)
      kPai(1) = kEndWij
      Do iSym = 2, nSym
         kPai(iSym) = kPai(iSym-1) + nVir(iSym-1)*nOcc(iSym-1)
         lPai = lPai + nVir(iSym)*nOcc(iSym)
      End Do
      kEndPai = kPai(1) + lPai
      Call FZero(Wrk(kPai(1)),lPai)

*     Allocate memory for Wai-Vector
*     ------------------------------

      lWai = nVir(1)*nOcc(1)
      kWai(1) = kEndPai
      Do iSym = 2, nSym
         kWai(iSym) = kWai(iSym-1) + nVir(iSym-1)*nOcc(iSym-1)
         lWai = lWai + nVir(iSym)*nOcc(iSym)
      End Do
      kEndWai = kWai(1) + lWai
      Call FZero(Wrk(kWai(1)),lWai)

*     Allocate memory for PaK-Vector
*     ------------------------------

      lPaK = nVir(1)*nFro(1)
      kPaK(1) = kEndWai
      Do iSym = 2, nSym
         kPaK(iSym) = kPaK(iSym-1) + nVir(iSym-1)*nFro(iSym-1)
         lPaK = lPaK + nVir(iSym)*nFro(iSym)
      End Do
      kEndPaK = kPaK(1) + lPaK
      Call FZero(Wrk(kPaK(1)),lPaK)

*     Allocate memory for WaK-Vector
*     ------------------------------

      lWaK = nVir(1)*nFro(1)
      kWaK(1) = kEndPaK
      Do iSym = 2, nSym
         kWaK(iSym) = kWaK(iSym-1) + nVir(iSym-1)*nOcc(iSym-1)
         lWaK = lWaK + nVir(iSym)*nOcc(iSym)
      End Do
      kEndWaK = kWaK(1) + lWaK
      Call FZero(Wrk(kWaK(1)),lWaK)


*     Allocate memory for PiK-vector (occ-fro)
*     ----------------------------------------

      lPiK = nOcc(1)*nFro(1)
      kPiK(1) = kEndWaK
      Do iSym = 2, nSym
         kPiK(iSym) = kPiK(iSym-1) + nOcc(iSym-1)*nFro(iSym-1)
         lPiK = lPiK + nOcc(iSym)*nFro(iSym)
      End Do
      kEndPiK = kPiK(1) + lPiK
      Call FZero(Wrk(kPiK(1)),lPiK)

*     Allocate memory for WiK-vector (occ-fro)
*     ----------------------------------------

      lWiK = nOcc(1)*nFro(1)
      kWiK(1) = kEndPiK
      Do iSym = 2, nSym
         kWiK(iSym) = kWiK(iSym-1) + nOcc(iSym-1)*nFro(iSym-1)
         lWiK = lWiK + nOcc(iSym)*nFro(iSym)
      End Do
      kEndWiK = kWiK(1) + lWiK
      Call FZero(Wrk(kWiK(1)),lWiK)

*     Allocate memory for WJK-vector (fro-fro)
*     ----------------------------------------

      lWJK = nFro(1)*nFro(1)
      kWJK(1) = kEndWiK
      Do iSym = 2, nSym
         kWJK(iSym) = kWJK(iSym-1) + nFro(iSym-1)*nFro(iSym-1)
         lWJK = lWJK + nFro(iSym)*nFro(iSym)
      End Do
      kEndWJK = kWJK(1) + lWJK
      Call FZero(Wrk(kWJK(1)),lWJK)

*     Allocate memory for Lagr-vector
*     ------------------------------

      lLagr = nOcc(1)*nVir(1)
      kLagr(1) = kEndWJK
      Do iSym = 2, nSym
         kLagr(iSym) = kLagr(iSym-1) + nOcc(iSym-1)*nVir(iSym-1)
         lLagr = lLagr + nOcc(iSym)*nVir(iSym)
      End Do
      kEndLagr = kLagr(1) + lLagr
      Call FZero(Wrk(kLagr(1)),lLagr)

*     Allocate memory for FrozenLagr-vector
*     -------------------------------------

      lFLagr = nFro(1)*nVir(1)
      kFLagr(1) = kEndLagr
      Do iSym = 2, nSym
         kFLagr(iSym) = kFLagr(iSym-1) + nFro(iSym-1)*nVir(iSym-1)
         lFLagr = lFLagr + nFro(iSym)*nVir(iSym)
      End Do
      kEndFLagr = kFLagr(1) + lFLagr
      Call FZero(Wrk(kFlagr(1)),lFLagr)

      nVec = Min(maxvalue,Max(nMP2VecMax,NumChoMax))
      If(nVec .lt. 1) Then
         Call ChoMP2_Quit(SecNam,'Insufficient memory','[1]')
      End If

*     Allocate memory for X^KJ-vector
*     -------------------------------

      lX = nVec*nMp2VecMax
      kX = kEndFLagr
      kEndX = kX + lX
*
      Do iSym = 1, nSym

*        Allocate memory for Ria-vectors
*        -------------------------------
         lRia = nMoMo(iSym,iVecOV)*nVec
         kRia = kEndX
         kEndRia = kRia + lRia
         kLia = kRia
*
         lRia2 = nMoMo(iSym,iVecOV)*nVec
         kRia2 = kEndRia
         kEndRia2 = kRia2 + lRia2

*        Allocate memory for L-vectors
*        -----------------------------
         lLKa = nMoMo(iSym,iVecFV)*nVec
         kLKa = kEndRia2
         kEndLKa = kLKa + lLKa

         lLab = nMoMo(iSym,iVecVV)*nVec
         kLab = kEndLKa
         kEndLab = kLab + lLab

         lLij = nMoMo(iSym,iVecOO)*nVec
         kLij = kEndLab
         kEndLij = kLij + lLij

         lLiK = nMoMo(iSym,iVecOF)*nVec
         kLiK = kEndLij
         kEndLiK = kLiK + lLiK

*        Allocate memory for U-vector
*        ----------------------------

         lU = nMoMo(iSym,iVecOV)*nVec
         kU = kEndLiK
         kEndU = kU + lU

*        Setup batch over amplitude vectors.
*        -----------------------------------
         nBatR = (nMP2Vec(iSym)-1)/nVec + 1
         nBatL = (NumCho(iSym)-1)/nVec + 1
         If((nMoMo(iSym,iVecOV) .gt. 0).and.(nMp2Vec(iSym).gt.0)) Then
*
*           Open the File for reordered amplitude vectors
*           ---------------------------------------------
            Call ChoMP2_OpenF(1,iTypR,iSym)
            Call ChoMP2_OpenF(1,iTypL,iSym)

            iSeed = 7
            LuUVec = IsFreeUnit(iSeed)
            Write(Fname,'(A4,I1)') 'TMPV',1
            Call DaName_MF_WA(LuUVec,Fname)
*
            iSeed = 7
            LuVVec = IsFreeUnit(iSeed)
            Write(Fname,'(A4,I1)') 'TMPV',2
            Call DaName_MF_WA(LuVVec,Fname)
*
*           Calculate Intermediate vectors U
*           --------------------------------
            Do kBat = 1, nBatR
               Call FZero(Wrk(kX),lX)
               If (kBat .eq. nBatR) Then
                  NumVecK = nMP2Vec(iSym) - nVec*(nBatR-1)
               Else
                  NumVecK = nVec
               End If
               kVec = nVec*(kBat-1) + 1

*              Read Amplitude vectors from kBat
*              --------------------------------
               iOpt = 2
               lTot = nMoMo(iSym,iVecOV)*NumVecK
               iAdr = nMoMo(iSym,iVecOV)*(kVec-1) + 1
               Call dDaFile(lUnit_F(iSym,iTypR),iOpt,Wrk(kRia),lTot,
     &                      iAdr)
*
               Do jBat = 1, nBatR
                  If (jBat .eq. nBatR) Then
                     NumVecJ = nMP2Vec(iSym) - nVec*(nBatR-1)
                  Else
                     NumVecJ = nVec
                  End If
                  jVec = nVec*(jBat-1) + 1
*
*                 Read Amplitude vectors from jBat
*                 --------------------------------
                  iOpt = 2
                  lTot = nMoMo(iSym,iVecOV)*NumVecJ
                  iAdr = nMoMo(iSym,iVecOV)*(jVec-1) + 1
                  Call dDaFile(lUnit_F(iSym,iTypR),iOpt,Wrk(kRia2),lTot,
     &                         iAdr)

*                 Construct X^JK-vector
*                 ---------------------
                 iOffX = NumVecK*(jVec-1)
                 Call dGemm_('T','N',NumVecK,NumVecJ,nMoMo(iSym,iVecOV),
     &                       1.0d0, Wrk(kRia),nMoMo(iSym,iVecOV),
     &                       Wrk(kRia2),nMoMo(iSym,iVecOV),0.0d0,
     &                       Wrk(kX+iOffX),NumVecK)

               End Do
*
               iOpt = 1
               lTot = nMP2Vec(iSym)*NumVecK
               iAdr = 1 + nMP2Vec(iSym)*(kVec-1)
               Call dDaFile(LuVVec,iOpt,Wrk(kX),lTot,iAdr)
            End Do
            Do kBat = 1, nBatR
               If (kBat .eq. nBatR) Then
                  NumVecK = nMP2Vec(iSym) - nVec*(nBatR-1)
               Else
                  NumVecK = nVec
               End If
               kVec = nVec*(kBat-1) + 1
*
               iOpt = 2
               lTot = nMP2Vec(iSym)*NumVecK
               iAdr = 1 + nMP2Vec(iSym)*(kVec-1)
               Call dDaFile(LuVVec,iOpt,Wrk(kX),lTot,iAdr)

               Do jBat = 1, nBatR
                  If (jBat .eq. nBatR) Then
                     NumVecJ = nMP2Vec(iSym) - nVec*(nBatR-1)
                  Else
                     NumVecJ = nVec
                  End If
                  jVec = nVec*(jBat-1) + 1

*                 Read Amplitude vectors from jBat
*                 --------------------------------
                  iOpt = 2
                  lTot = nMoMo(iSym,iVecOV)*NumVecJ
                  iAdr = nMoMo(iSym,iVecOV)*(jVec-1) + 1
                  Call dDaFile(lUnit_F(iSym,iTypR),iOpt,Wrk(kRia),lTot,
     &                         iAdr)
*
                 iOffX = NumVecK*(jVec-1)
                 Fac = X(min((jBat-1),1))
                 Call dGemm_('N','T',nMoMo(iSym,iVecOV),NumVecK,NumVecJ,
     &                       1.0d0,Wrk(kRia),nMoMo(iSym,iVecOV),
     &                       Wrk(kX+iOffX),NumVecK,Fac,
     &                       Wrk(kU),nMoMo(iSym,iVecOV))
               End Do
*
               iOpt = 1
               lTot = nMoMo(iSym,iVecOV)*NumVecK
               iAdr = 1 + nMoMo(iSym,iVecOV)*(kVec-1)
               Call dDaFile(LuUVec,iOpt,Wrk(kU),
     &                      lTot,iAdr)
            End Do
*
*           Calculate "Coulomb"-contributions to Densities from U
*           -----------------------------------------------------
            Do kBat = 1,nBatR
               If (kBat .eq. nBatR) Then
                  NumVecK = nMP2Vec(iSym) - nVec*(nBatR-1)
               Else
                  NumVecK = nVec
               End If
               kVec = nVec*(kBat-1) + 1

*              Read Amplitude vectors from kBat
*              --------------------------------
               iOpt = 2
               lTot = nMoMo(iSym,iVecOV)*NumVecK
               iAdr = nMoMo(iSym,iVecOV)*(kVec-1) + 1
               Call dDaFile(lUnit_F(iSym,iTypR),iOpt,Wrk(kRia),lTot,
     &                      iAdr)

*              Load intermediate vectors U^K_ib
*              --------------------------------
               iOpt = 2
               lTot = nMoMo(iSym,iVecOV)*NumVecK
               iAdr = nMoMo(iSym,iVecOv)*(kVec-1) + 1
               Call dDaFile(LuUVec,iOpt,Wrk(kU),lTot,iAdr)

*              Calculate the "Coulomb" Contribution to Pab
*              -------------------------------------------
               iOff1 = 0
               Do iSymI = 1, nSym
                  iSymA = MulD2h(iSym,iSymI)
                  If(nOcc(iSymI)*nVir(iSymA).eq.0) Go To 121
                  Do iI = 1, nOcc(iSymI)
                     iOff = (iI-1)*nVir(iSymA) + iOff1
                     Call dGemm_('N','T',nVir(iSymA),nVir(iSymA),
     &                          NumVecK, 4.0d0,
     &                          Wrk(kRia+iOff),nMoMo(iSym,iVecOV),
     &                          Wrk(kU+iOff),nMoMo(iSym,iVecOV),1.0d0,
     &                          Wrk(kPab(iSymA)),nVir(iSymA))
                  End Do
 121              Continue
                  iOff1 = iOff1 + nOcc(iSymI)*nVir(iSymA)
               End Do
*              Calculate the "Coulomb" Contribution to Pik
*              -------------------------------------------
*
               Do kVec1 = 1, NumVecK
                  iOff1 = 0
                  Do iSymK = 1, nSym
                     iSymI = iSymK
                     iSymB = MulD2h(iSymK,iSym)
                     If(nOcc(iSymK)*nVir(iSymB) .eq. 0) Go To 122
                     iOff = nMoMo(iSym,iVecOV)*(kVec1-1) + iOff1
                     Call dGemm_('T','N', nOcc(iSymI), nOcc(iSymK),
     &                          nVir(iSymB), -4.0d0,
     &                          Wrk(kU+iOff),nVir(iSymB),
     &                          Wrk(kRia+iOff),nVir(iSymB),1.0d0,
     &                          Wrk(kPij(iSymK)),nOcc(iSymI))
 122                 Continue
                     iOff1 = iOff1 + nOcc(iSymI)*nVir(iSymB)
                  End Do
               End Do

            End Do
*

*           Calculate Intermediate vectors U*
*           ---------------------------------
            Do kBat = 1, nBatL
               Call FZero(Wrk(kX),lX)
               If (kBat .eq. nBatL) Then
                  NumVecK = NumCho(iSym) - nVec*(nBatL-1)
               Else
                  NumVecK = nVec
               End If
               kVec = nVec*(kBat-1) + 1

*              Read Integral vectors from kBat
*              --------------------------------
               iOpt = 2
               lTot = nMoMo(iSym,iVecOV)*NumVecK
               iAdr = 1 + nMoMo(iSym,iVecOV)*(kVec-1) +
     &                iAdrOff(iSym,iVecOV)
               Call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLia),lTot,
     &                      iAdr)
*
               Do jBat = 1, nBatR
                  If (jBat .eq. nBatR) Then
                     NumVecJ = nMP2Vec(iSym) - nVec*(nBatR-1)
                  Else
                     NumVecJ = nVec
                  End If
                  jVec = nVec*(jBat-1) + 1
*
*                 Read Amplitude vectors from jBat
*                 --------------------------------
                  iOpt = 2
                  lTot = nMoMo(iSym,iVecOV)*NumVecJ
                  iAdr = nMoMo(iSym,iVecOV)*(jVec-1) + 1
                  Call dDaFile(lUnit_F(iSym,iTypR),iOpt,Wrk(kRia2),lTot,
     &                         iAdr)

*                 Construct X^JK-vector
*                 ---------------------
                 iOffX = NumVecK*(jVec-1)
                 Call dGemm_('T','N',NumVecK,NumVecJ,nMoMo(iSym,iVecOV),
     &                       1.0d0, Wrk(kLia),nMoMo(iSym,iVecOV),
     &                       Wrk(kRia2),nMoMo(iSym,iVecOV),0.0d0,
     &                       Wrk(kX+iOffX),NumVecK)
               End Do
*
               iOpt = 1
               lTot = nMP2Vec(iSym)*NumVecK
               iAdr = 1 + nMP2Vec(iSym)*(kVec-1)
               Call dDaFile(LuVVec,iOpt,Wrk(kX),lTot,iAdr)
            End Do
            Do kBat = 1, nBatL
               If (kBat .eq. nBatL) Then
                  NumVecK = NumCho(iSym) - nVec*(nBatL-1)
               Else
                  NumVecK = nVec
               End If
               kVec = nVec*(kBat-1) + 1
*
               iOpt = 2
               lTot = nMP2Vec(iSym)*NumVecK
               iAdr = 1 + nMP2Vec(iSym)*(kVec-1)
               Call dDaFile(LuVVec,iOpt,Wrk(kX),lTot,iAdr)
*
               Do jBat = 1, nBatR
                  If (jBat .eq. nBatR) Then
                     NumVecJ = nMP2Vec(iSym) - nVec*(nBatR-1)
                  Else
                     NumVecJ = nVec
                  End If
                  jVec = nVec*(jBat-1) + 1

*                 Read Amplitude vectors from jBat
*                 --------------------------------
                  iOpt = 2
                  lTot = nMoMo(iSym,iVecOV)*NumVecJ
                  iAdr = nMoMo(iSym,iVecOV)*(jVec-1) + 1
                  Call dDaFile(lUnit_F(iSym,iTypR),iOpt,Wrk(kRia),lTot,
     &                         iAdr)
*
                 iOffX = NumVecK*(jVec-1)
                 Fac = X(min((jBat-1),1))
                 Call dGemm_('N','T',nMoMo(iSym,iVecOV),NumVecK,NumVecJ,
     &                       1.0d0,Wrk(kRia),nMoMo(iSym,iVecOV),
     &                       Wrk(kX+iOffX),NumVecK,Fac,
     &                       Wrk(kU),nMoMo(iSym,iVecOV))
               End Do
*
               iOpt = 1
               lTot = nMoMo(iSym,iVecOV)*NumVecK
               iAdr = 1 + nMoMo(iSym,iVecOV)*(kVec-1)
               Call dDaFile(LuUVec,iOpt,Wrk(kU),lTot,iAdr)
            End Do

*           Calculate "Coulomb"-contributions to Densities and
*           the Lagrangian from U*
*           -----------------------------------------------------
            Do kBat = 1,nBatL
*
               If (kBat .eq. nBatL) Then
                  NumVecK = NumCho(iSym) - nVec*(nBatL-1)
               Else
                  NumVecK = nVec
               End If
               kVec = nVec*(kBat-1) + 1

*              Read Integral vectors(FV) from kBat
*              --------------------------------
               iOpt = 2
               lTot = nMoMo(iSym,iVecFV)*NumVecK
               iAdr = 1 + nMoMo(iSym,iVecFV)*(kVec-1) +
     &                 iAdrOff(iSym,iVecFV)
               Call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLKa),lTot,
     &                      iAdr)
*
*              Read Integral vectors(VV) from kBat
*              --------------------------------
               iOpt = 2
               lTot = nMoMo(iSym,iVecVV)*NumVecK
               iAdr = 1 + nMoMo(iSym,iVecVV)*(kVec-1) +
     &                 iAdrOff(iSym,iVecVV)
               Call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLab),lTot,
     &                      iAdr)

*              Read Integral vectors(OO) from kBat
*              --------------------------------
               iOpt = 2
               lTot = nMoMo(iSym,iVecOO)*NumVecK
               iAdr = 1 + nMoMo(iSym,iVecOO)*(kVec-1) +
     &                 iAdrOff(iSym,iVecOO)
               Call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLij),lTot,
     &                      iAdr)

*              Read Integral vectors(OF) from kBat
*              --------------------------------
               iOpt = 2
               lTot = nMoMo(iSym,iVecOF)*NumVecK
               iAdr = 1 + nMoMo(iSym,iVecOF)*(kVec-1) +
     &                 iAdrOff(iSym,iVecOF)
               Call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLiK),lTot,
     &                      iAdr)

*              Read Integral vectors(OV) from kBat
*              --------------------------------
               iOpt = 2
               lTot = nMoMo(iSym,iVecOV)*NumVecK
               iAdr = 1 + nMoMo(iSym,iVecOV)*(kVec-1) +
     &                 iAdrOff(iSym,iVecOV)
               Call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLia),lTot,
     &                      iAdr)

*              Load intermediate vectors
*              -------------------------
               iOpt = 2
               lTot = nMoMo(iSym,iVecOV)*NumVecK
               iAdr = nMoMo(iSym,iVecOv)*(kVec-1) + 1
               Call dDaFile(LuUVec,iOpt,Wrk(kU),lTot,iAdr)

*        Calculate the "Coulomb" Contribution to Froz-occupied density PiK
*        -----------------------------------------------------------------
*
               Do kVec1 = 1, NumVecK
                  iOffU1 = 0
                  iOffL1 = 0
                  Do iSymK = 1, nSym
                     iSymI = iSymK
                     iSymB = MulD2h(iSymK,iSym)
                     If(nOcc(iSymI)*nFro(iSymK)*nVir(iSymB) .eq. 0)
     &                      Go To 131
                    iOffU = nMoMo(iSym,iVecOV)*(kVec1-1) + iOffU1
                    iOffL = nMoMo(iSym,iVecFV)*(kVec1-1) + iOffL1
                    Call dGemm_('T','N', nOcc(iSymI), nFro(iSymK),
     &                          nVir(iSymB), -4.0d0,
     &                          Wrk(kU+iOffU),nVir(iSymB),
     &                          Wrk(kLKa+iOffL),nVir(iSymB),1.0d0,
     &                          Wrk(kPiK(iSymK)),nOcc(iSymI))
 131                 Continue
                     iOffU1 = iOffU1 + nOcc(iSymI)*nVir(iSymB)
                     iOffL1 = iOffL1 + nFro(iSymK)*nVir(iSymB)
                  End Do
               End Do

*              Calculate the "Coulomb" contribution to Wij (froz-occ)
*              ------------------------------------------------------

               Do kVec1 = 1, NumVecK
                  iOffU1 = 0
                  iOffL1 = 0
                  Do iSymK = 1, nSym
                     iSymI = iSymK
                     iSymB = MulD2h(iSymK,iSym)
                     If(nOcc(iSymI)*nFro(iSymK)*nVir(iSymB) .eq. 0)
     &                      Go To 135
                    iOffU = nMoMo(iSym,iVecOV)*(kVec1-1) + iOffU1
                    iOffL = nMoMo(iSym,iVecFV)*(kVec1-1) + iOffL1
                    Call dGemm_('T','N', nOcc(iSymI), nFro(iSymK),
     &                          nVir(iSymB), 4.0d0,
     &                          Wrk(kU+iOffU),nVir(iSymB),
     &                          Wrk(kLKa+iOffL),nVir(iSymB),1.0d0,
     &                          Wrk(kWiK(iSymK)),nOcc(iSymI))
 135                 Continue
                     iOffU1 = iOffU1 + nOcc(iSymI)*nVir(iSymB)
                     iOffL1 = iOffL1 + nFro(iSymK)*nVir(iSymB)
                  End Do
               End Do


*              Calculate "Coulomb"-contribution to Lagr(IV)
*              --------------------------------------------
               Do kVec1 = 1, NumVecK
                  iOffU1 = 0
                  iOffL1 = 0
                  Do iSymA = 1, nSym
                     iSymI = iSymA
                     iSymB = MulD2h(iSymA,iSym)
                     If(nOcc(iSymI)*nVir(iSymA)*nVir(iSymB) .eq. 0)
     &                      Go To 132
                     iOffU = nMoMo(iSym,iVecOV)*(kVec1-1) + iOffU1
                     iOffL = nMoMo(iSym,iVecVV)*(kVec1-1) + iOffL1
                     Call dGemm_('T','N', nVir(iSymA), nOcc(iSymI),
     &                          nVir(iSymB), 4.0d0,
     &                          Wrk(kLab+iOffL),nVir(iSymB),
     &                          Wrk(kU+iOffU),nVir(iSymB),1.0d0,
     &                          Wrk(kLagr(iSymA)),nVir(iSymA))
*
 132                 Continue
                     iOffU1 = iOffU1 + nOcc(iSymI)*nVir(iSymB)
                     iOffL1 = iOffL1 + nVir(iSymA)*nVir(iSymB)
                  End Do
               End Do

*              Calculate "Coulomb"-contribution to Lagr(III)
*              ---------------------------------------------
               iOffU1 = 0
               iOffL1 = 0
               Do iSymJ = 1, nSym
                  iSymA = MulD2h(iSym,iSymJ)
                  iSymI = iSymA
                  If(nOcc(iSymI)*nVir(iSymA)*nOcc(iSymJ).eq.0) Go To 133
                  Do iJ = 1, nOcc(iSymJ)
                     iOffU = (iJ-1)*nVir(iSymA) + iOffU1
                     iOffL = (iJ-1)*nOcc(iSymI) + iOffL1
                     Call dGemm_('N','T',nVir(iSymA),nOcc(iSymI),
     &                         NumVecK, -4.0d0,
     &                         Wrk(kU+iOffU),nMoMo(iSym,iVecOV),
     &                         Wrk(kLij+iOffL),nMoMo(iSym,iVecOO),1.0d0,
     &                         Wrk(kLagr(iSymA)),nVir(iSymA))
                  End Do
 133              Continue
                  iOffU1 = iOffU1 + nOcc(iSymJ)*nVir(iSymA)
                  iOffL1 = iOffL1 + nOcc(iSymJ)*nOcc(iSymI)
               End Do

*              Calculate "Coulomb"-contribution to Wai
*              ---------------------------------------
               iOffU1 = 0
               iOffL1 = 0
               Do iSymJ = 1, nSym
                  iSymA = MulD2h(iSym,iSymJ)
                  iSymI = iSymA
                  If(nOcc(iSymI)*nVir(iSymA)*nOcc(iSymJ).eq.0) GoTo 1331
                  Do iJ = 1, nOcc(iSymJ)
                     iOffU = (iJ-1)*nVir(iSymA) + iOffU1
                     iOffL = (iJ-1)*nOcc(iSymI) + iOffL1
                     Call dGemm_('N','T',nVir(iSymA),nOcc(iSymI),
     &                         NumVecK, 8.0d0,
     &                         Wrk(kU+iOffU),nMoMo(iSym,iVecOV),
     &                         Wrk(kLij+iOffL),nMoMo(iSym,iVecOO),1.0d0,
     &                         Wrk(kWai(iSymA)),nVir(iSymA))
                  End Do
 1331             Continue
                  iOffU1 = iOffU1 + nOcc(iSymJ)*nVir(iSymA)
                  iOffL1 = iOffL1 + nOcc(iSymJ)*nOcc(iSymI)
               End Do


*              Calculate "Coulomb"-contribution to Wab
*              ---------------------------------------
               iOff1 = 0
               Do iSymI = 1, nSym
                  iSymA = MulD2h(iSym,iSymI)
                  If(nOcc(iSymI)*nVir(iSymA).eq.0) Go To 151
                  Do iI = 1, nOcc(iSymI)
                     iOff = (iI-1)*nVir(iSymA)+iOff1
                     Call dGemm_('N','T',nVir(iSymA),nVir(iSymA),
     &                          NumVecK, 4.0d0,
     &                          Wrk(kLia+iOff),nMoMo(iSym,iVecOV),
     &                          Wrk(kU+iOff),nMoMo(iSym,iVecOV),1.0d0,
     &                          Wrk(kWab(iSymA)),nVir(iSymA))
                  End Do
 151              Continue
                  iOff1 = iOff1 + nOcc(iSymI)*nVir(iSymA)
               End Do

*              Calculate "Coulomb"-contribution to Wij (Wik)
*              ---------------------------------------------
               Do kVec1 = 1, NumVecK
                  iOff1 = 0
                  Do iSymK = 1, nSym
                     iSymI = iSymK
                     iSymB = MulD2h(iSymK,iSym)
                     If(nOcc(iSymK)*nVir(iSymB) .eq. 0) Go To 152
                     iOff = nMoMo(iSym,iVecOV)*(kVec1-1) + iOff1
                     Call dGemm_('T','N', nOcc(iSymI), nOcc(iSymK),
     &                          nVir(iSymB), 4.0d0,
     &                          Wrk(kU+iOff),nVir(iSymB),
     &                          Wrk(kLia+iOff),nVir(iSymB),1.0d0,
     &                          Wrk(kWij(iSymK)),nOcc(iSymI))
 152                 Continue
                  End Do
               End Do


*     Calculate "Coulomb"-contribution to the Frozen Lagr(III)
*     --------------------------------------------------------
               iOffU1 = 0
               iOffL1 = 0
               Do iSymJ = 1, nSym
                  iSymA = MulD2h(iSym,iSymJ)
                  iSymI = iSymA
                  If(nFro(iSymI)*nVir(iSymA)*nOcc(iSymJ).eq.0) Go To 134
                  Do iJ = 1, nOcc(iSymJ)
                     iOffU = (iJ-1)*nVir(iSymA) + iOffU1
                     iOffL = (iJ-1)*nFro(iSymI) + iOffL1
                     Call dGemm_('N','T',nVir(iSymA),nFro(iSymI),
     &                         NumVecK, -4.0d0,
     &                         Wrk(kU+iOffU),nMoMo(iSym,iVecOV),
     &                         Wrk(kLiK+iOffL),nMoMo(iSym,iVecOF),1.0d0,
     &                         Wrk(kFLagr(iSymA)),nVir(iSymA))
                  End Do
 134              Continue
                  iOffU1 = iOffU1 + nOcc(iSymJ)*nVir(iSymA)
                  iOffL1 = iOffL1 + nOcc(iSymJ)*nFro(iSymI)
               End Do


*     Calculate "Coulomb"-contribution to the WaK (vir-froz)
*     --------------------------------------------------------
               iOffU1 = 0
               iOffL1 = 0
               Do iSymJ = 1, nSym
                  iSymA = MulD2h(iSym,iSymJ)
                  iSymI = iSymA
                  If(nFro(iSymI)*nVir(iSymA)*nOcc(iSymJ).eq.0) GoTo 1341
                  Do iJ = 1, nOcc(iSymJ)
                     iOffU = (iJ-1)*nVir(iSymA) + iOffU1
                     iOffL = (iJ-1)*nFro(iSymI) + iOffL1
                     Call dGemm_('N','T',nVir(iSymA),nFro(iSymI),
     &                         NumVecK, 8.0d0,
     &                         Wrk(kU+iOffU),nMoMo(iSym,iVecOV),
     &                         Wrk(kLiK+iOffL),nMoMo(iSym,iVecOF),1.0d0,
     &                         Wrk(kWaK(iSymA)),nVir(iSymA))
                  End Do
 1341             Continue
                  iOffU1 = iOffU1 + nOcc(iSymJ)*nVir(iSymA)
                  iOffL1 = iOffL1 + nOcc(iSymJ)*nFro(iSymI)
               End Do
            End Do

            Call ChoMP2_OpenF(iClos,iTypR,iSym)
            Call ChoMP2_OpenF(iClos,iTypL,iSym)
            Call DaClos(LuUVec)
            Call DaClos(LuVVec)

         End If                 ! vector check
      End Do !iSym

***********************************************************
*     --------------------------------------------------  *
*     Calculate the exchange part needing U_ic^K-vectors  *
*     --------------------------------------------------  *
***********************************************************
*     Calculate the max number of orbitals of a specific type
*     -------------------------------------------------------
      nFroMax = 0
      nVirMax = 0
      nOccMax = 0
      Do iSym = 1, nSym
         nFroMax = Max(nFroMax,nFro(iSym))
         nVirMax = Max(nVirMax,nVir(iSym))
         nOccMax = Max(nOccMax,nOcc(iSym))
      End Do



      If(nVirMax .gt. 1) Then
         nB = Min(nVirMax,maxvalue)
      Else
         nB = nVirMax
      End If

*     Allocate memory for T_[i]j^[b]c
*     -------------------------------
      lAmp = nVirMax*nB*nOccMax
      kAmp = kEndFLagr
      kEndAmp = kAmp + lAmp
      Call FZero(Wrk(kAmp),lAmp)

*
*     Allocate memory for Rjc^K
*     -------------------------
      lRjc = nVirMax*nVec
      kRjc = kEndAmp
      kEndRjc = kRjc+lRjc

*     Allocate memory for Rib^K
*     -------------------------
      lRib = nOccMax*nVec
      kRib = kEndRjc
      kEndRib = kRib + lRib


*     Open files for Reordered R-vectors
*     ----------------------------------
      iSeed = 7
      LuRInv(1) = IsFreeUnit(iSeed)
      Write(Fname,'(A4,I1)') 'TMPV',4
      Call DaName_MF_WA(LuRInv(1),Fname)
*
      iSeed = 8
      LuRInv(2) = IsFreeUnit(iSeed)
      Write(Fname,'(A4,I1)') 'TMPV',5
      Call DaName_MF_WA(LuRInv(2),Fname)

*
      Do jSym = 1, nSym
         If((nMoMo(jSym,iVecOV) .le. 0) .or.
     &      (nMp2Vec(jSym)*NumCho(jSym)).le.0)
     &       Go To 10

*        Allocate memory for Ric^J
*        -------------------------
         lRic = nMoMo(jSym,iVecOV)*nVec
         kRic = kEndRib
         kEndRic = kRic + lRic

*        Allocate memory for Lic^J
*        -------------------------
         lLic = nMoMo(jSym,iVecOV)*nVec
         kLic = kEndRic
         kEndLic = kLic + lLic

*        Allocate memory for Lab^J
*        -------------------------
         lLab = nMoMo(jSym,iVecVV)*nVec
         kLab = kEndLic
         kEndLab = kLab + lLab

*        Allocate memory for Lji^J
*        -------------------------
         lLji = nMoMo(jSym,iVecOO)*nVec
         kLji = kEndLab
         kEndLji = kLji + lLji

*        Allocate memory for LKa^J
*        -------------------------
         lLKa = nMoMo(jSym,iVecFV)*nVec
         kLKa = kEndLji
         kEndLKa = kLKa + lLKa

*        Allocate memory for LiK^J
*        -------------------------
         lLiK = nMoMo(jSym,iVecOF)*nVec
         kLiK = kEndLKa
         kEndLiK = kLiK + lLiK

*        Allocate memory for Ujb^J
*        -------------------------
         lU = nMoMo(jSym,iVecOV)*nVec
         kU = kEndLiK
         kEndU = kU+lU

*        Allocate memory for Vjb^J
*        -------------------------
         lV = nMoMo(jSym,iVecOV)*nVec
         kV = kEndU

         nBatR = (nMP2Vec(jSym)-1)/nVec + 1
         nBatL = (NumCho(jSym)-1)/nVec + 1
*
*        Open the File for reordered amplitude vectors
*        ---------------------------------------------
         Call ChoMP2_OpenF(1,iTypR,jSym)
         Call ChoMP2_OpenF(1,iTypL,jSym)

         iSeed = 7
         LuUVec = IsFreeUnit(iSeed)
         Write(Fname,'(A4,I1)') 'TMPV',1
         Call DaName_MF_WA(LuUVec,Fname)

         iSeed = 7
         LuVVec = IsFreeUnit(iSeed)
         Write(Fname,'(A4,I1)') 'TMPV',2
         Call DaName_MF_WA(LuVVec,Fname)

         nBatMax = Max(nBatR,nBatL)
         Do jBat = 1, nBatMax
            Call fZero(Wrk(kU),lU)
            Call fZero(Wrk(kV),lV)
            If (jBat .eq. nBatR) Then
               NumRVecJ = nMP2Vec(jSym) - nVec*(nBatR-1)
            Else If(jBat .gt. nBatR) Then
               NumRVecJ = 0
            Else
               NumRVecJ = nVec
            End If
            If (jBat .eq. nBatL) Then
               NumLVecJ = NumCho(jSym) - nVec*(nBatL-1)
            Else If(jBat .gt. nBatL) Then
               NumLVecJ = 0
            Else
               NumLVecJ = nVec
            End If
            jVec = nVec*(jBat-1) + 1

*           Read Amplitude vectors from kBat
*           --------------------------------
            If(NumRVecJ .gt. 0) Then
               iOpt = 2
               lTot = nMoMo(jSym,iVecOV)*NumRVecJ
               iAdr = nMoMo(jSym,iVecOV)*(jVec-1) + 1
               Call dDaFile(lUnit_F(jSym,iTypR),iOpt,Wrk(kRic),lTot,
     &                      iAdr)
            End If
            If(NumLVecJ .gt. 0) Then
               iOpt = 2
               lTot = nMoMo(jSym,iVecOV)*NumLVecJ
               iAdr = nMoMo(jSym,iVecOV)*(jVec-1) + 1
               Call dDaFile(lUnit_F(jSym,iTypL),iOpt,Wrk(kLic),lTot,
     &                      iAdr)
            End If

            Do iSymJ = 1, nSym

               iSymB = MulD2h(jSym,iSymJ)
               If(nVir(iSymB) .eq. 0) Go To 202
               nBBlock = (nVir(iSymB)-1)/nB + 1

               Do iSymI = 1, nSym
                  iSymC = MulD2h(jSym,iSymI)
                  iSymJC = MulD2h(iSymJ,iSymC)
                  NumIC = nOcc(iSymI)*nVir(iSymC)
                  If(nOcc(iSymI)*nVir(iSymC)*(NumRVecJ+NumLVecJ) .eq. 0)
     &              Go To 201
                  Do iJ = 1, nOcc(iSymJ)
                     Do iBBlock = 1, nBBlock
                        If(iBBlock.eq.nBBlock) Then
                           NumB = nVir(iSymB) - nB*(nBBlock-1)
                        Else
                           NumB = nB
                        End If
                        iB1 = nB*(iBBlock-1) + 1

                        Do kBat = 1, nBatR
                           If (kBat .eq. nBatR) Then
                              NumVecK = nMP2Vec(iSymJC) - nVec*(nBatR-1)
                           Else
                              NumVecK = nVec
                           End If
                           kVec = nVec*(kBat-1) + 1

                           iOpt = 2
                           lTot = NumVecK*nVir(iSymC)
                           iAdr = iWork(ipAdrR1+iAdrVec(iSymC,iSymJ,iJ))
     &                          + (kVec-1)*nVir(iSymC)
                           Call dDaFile(LuRInv(1),iOpt,Wrk(kRjc),
     &                                  lTot,iAdr)
*
                           Do iBRel = 1, NumB
                              iB = iBrel + (iB1-1)
                              iOpt = 2
                              lTot = NumVecK*nOcc(iSymI)
                              iAdr = iWork(ipAdrR2+
     &                                     iAdrVec(iSymB,iSymI,iB))
     &                             + (kVec-1)*nOcc(iSymI)
*
                              Call dDaFile(LuRInv(2),iOpt, Wrk(kRib),
     &                                     lTot,iAdr)
                              iOffAmp =(iBrel-1)*nVir(iSymC)*nOcc(iSymI)
                              Fac = X(min((kBat-1),1))
                              Call dGemm_('N','T',nVir(iSymC),
     &                             nOcc(iSymI),
     &                             NumVecK, 1.0d0,Wrk(kRjc),
     &                             nVir(iSymC),
     &                             Wrk(kRib),nOcc(iSymI),Fac,
     &                             Wrk(kAmp+iOffAmp),nVir(iSymC))
                           End Do
                        End Do

*                       Calculate the contribution to U^J_jb
*                       ------------------------------------
                        If(NumRVecJ .gt. 0) Then
                           iOffRic = iT1am(iSymC,iSymI)
                           iOffU = iT1am(iSymB,iSymJ)+(iJ-1)*nVir(iSymB)
     &                           + (iB1-1)
                           Call dGemm_('T','N', NumB, NumRVecJ, NumIC,
     &                           1.0d0,Wrk(kAmp),NumIC,
     &                           Wrk(kRic+iOffRic),nMoMo(jSym,iVecOV),
     &                           1.0d0,Wrk(kU+iOffU),nMoMo(jSym,iVecOV))
                        End If
*                       Calculate the contribution to V^J_jb
*                       ------------------------------------
                        If(NumLVecJ .gt. 0) Then
                           iOffLic = iT1am(iSymC,iSymI)
                           iOffV = iT1am(iSymB,iSymJ)+(iJ-1)*nVir(iSymB)
     &                           + (iB1-1)
                           Call dGemm_('T','N', NumB, NumLVecJ, NumIC,
     &                           1.0d0,Wrk(kAmp),NumIC,
     &                           Wrk(kLic+iOffLic),nMoMo(jSym,iVecOV),
     &                           1.0d0,Wrk(kV+iOffV),nMoMo(jSym,iVecOV))
                        End If

                     End Do

                  End Do        !iJ
 201              Continue
               End Do           !iSymI
 202           Continue
            End Do              !iSymJ
*
            iOpt = 1
            lTot = nMoMo(jSym,iVecOV)*NumRVecJ
            iAdr = 1 + nMoMo(jSym,iVecOV)*(jVec-1)
            Call dDaFile(LuUVec,iOpt,Wrk(kU),
     &                   lTot,iAdr)
*
            iOpt = 1
            lTot = nMoMo(jSym,iVecOV)*NumLVecJ
            iAdr = 1 + nMoMo(jSym,iVecOV)*(jVec-1)
            Call dDaFile(LuVVec,iOpt,Wrk(kV),
     &                   lTot,iAdr)
         End Do                 !jBat
*
         Do jBat = 1, nBatMax
            If (jBat .eq. nBatR) Then
               NumRVecJ = nMP2Vec(jSym) - nVec*(nBatR-1)
            Else If(jBat .gt. nBatR) Then
               NumRVecJ = 0
            Else
               NumRVecJ = nVec
            End If
            If (jBat .eq. nBatL) Then
               NumLVecJ = NumCho(jSym) - nVec*(nBatL-1)
            Else If(jBat .gt. nBatL) Then
               NumLVecJ = 0
            Else
               NumLVecJ = nVec
            End If
            jVec = nVec*(jBat-1) + 1

*           Read Amplitude Rja
*           --------------------------------
            If(NumRVecJ .gt. 0) Then
               iOpt = 2
               lTot = nMoMo(jSym,iVecOV)*NumRVecJ
               iAdr = nMoMo(jSym,iVecOV)*(jVec-1) + 1
               Call dDaFile(lUnit_F(jSym,iTypR),iOpt,Wrk(kRic),lTot,
     &              iAdr)
            End If
*           Read Integrals R_ab^J from disk
*           -------------------------------
            If(NumLVecJ .gt. 0) Then
               iOpt = 2
               lTot = nMoMo(jSym,iVecVV)*NumLVecJ
               iAdr = 1 + nMoMo(jSym,iVecVV)*(jVec-1) +
     &                 iAdrOff(jSym,iVecVV)
               Call dDaFile(lUnit_F(jSym,iTypL),iOpt,Wrk(kLab),lTot,
     &              iAdr)
            End If

*           Read Integrals L_ab^J from disk
*           -------------------------------
            If(NumLVecJ .gt. 0) Then
               iOpt = 2
               lTot = nMoMo(jSym,iVecVV)*NumLVecJ
               iAdr = 1 + nMoMo(jSym,iVecVV)*(jVec-1) +
     &                 iAdrOff(jSym,iVecVV)
               Call dDaFile(lUnit_F(jSym,iTypL),iOpt,Wrk(kLab),lTot,
     &              iAdr)
            End If

*           Read L_ji^J-vectors from disk
*           -----------------------------
            If(NumLVecJ .gt. 0) Then
               iOpt = 2
               lTot = nMoMo(jSym,iVecOO)*NumLVecJ
               iAdr = 1 + nMoMo(jSym,iVecOO)*(jVec-1) +
     &                 iAdrOff(jSym,iVecOO)
               Call dDaFile(lUnit_F(jSym,iTypL),iOpt,Wrk(kLji),lTot,
     &              iAdr)
            End If

*           Read L_ic^J-vectors from disk
*           -----------------------------
            If(NumLVecJ .gt. 0) Then
               iOpt = 2
               lTot = nMoMo(jSym,iVecOV)*NumLVecJ
               iAdr = nMoMo(jSym,iVecOV)*(jVec-1) + 1
               Call dDaFile(lUnit_F(jSym,iTypL),iOpt,Wrk(kLic),lTot,
     &                      iAdr)
            End If

*           Read L_Ka^J-vectors from disk
*           -----------------------------
            If(NumLVecJ .gt. 0) Then
               iOpt = 2
               lTot = nMoMo(jSym,iVecFV)*NumLVecJ
               iAdr = 1 + nMoMo(jSym,iVecFV)*(jVec-1) +
     &                 iAdrOff(jSym,iVecFV)
               Call dDaFile(lUnit_F(jSym,iTypL),iOpt,Wrk(kLKa),lTot,
     &              iAdr)
            End If

*           Read L_Ka^J-vectors from disk
*           -----------------------------
            If(NumLVecJ .gt. 0) Then
               iOpt = 2
               lTot = nMoMo(jSym,iVecOF)*NumLVecJ
               iAdr = 1 + nMoMo(jSym,iVecOF)*(jVec-1) +
     &                 iAdrOff(jSym,iVecOF)
               Call dDaFile(lUnit_F(jSym,iTypL),iOpt,Wrk(kLiK),lTot,
     &              iAdr)
            End If

*           Read U_jb^J-vectors from disk
*           -----------------------------
            If(NumRVecJ .gt. 0) Then
               iOpt = 2
               lTot = nMoMo(jSym,iVecOV)*NumRVecJ
               iAdr = nMoMo(jSym,iVecOv)*(jVec-1) + 1
               Call dDaFile(LuUVec,iOpt,Wrk(kU),lTot,iAdr)
            End If


*           Read V_jb^J-vectors from disk
*           -----------------------------
            If(NumLVecJ .gt. 0) Then
               iOpt = 2
               lTot = nMoMo(jSym,iVecOV)*NumLVecJ
               iAdr = nMoMo(jSym,iVecOv)*(jVec-1) + 1
               Call dDaFile(LuVVec,iOpt,Wrk(kV),lTot,iAdr)
            End If


*           Calculate the "Exchange" Contribution to Pab
*           -------------------------------------------

            iOff1 = 0
            Do iSymJ = 1, nSym
               iSymA = MulD2h(jSym,iSymJ)
               If(nOcc(iSymJ)*nVir(iSymA)*NumRVecJ.eq.0) Go To 221
               Do iI = 1, nOcc(iSymJ)
                  iOff = (iI-1)*nVir(iSymA) + iOff1
                  Call dGemm_('N','T',nVir(iSymA),nVir(iSymA),
     &                       NumRVecJ, -2.0d0,
     &                       Wrk(kRic+iOff),nMoMo(jSym,iVecOV),
     &                       Wrk(kU+iOff),nMoMo(jSym,iVecOV),1.0d0,
     &                       Wrk(kPab(iSymA)),nVir(iSymA))
               End Do
 221           Continue
               iOff1 = iOff1 + nOcc(iSymJ)*nVir(iSymA)
            End Do

*           Calculate the "Exchange" Contribution to Pij
*           -------------------------------------------
*
            Do jVec1 = 1, NumRVecJ
               iOff1 = 0
               Do iSymJ = 1, nSym
                  iSymI = iSymJ
                  iSymB = MulD2h(iSymJ,jSym)
                  If(nOcc(iSymJ)*nVir(iSymB) .eq. 0) Go To 222
                  iOff = nMoMo(jSym,iVecOV)*(jVec1-1) + iOff1
                  Call dGemm_('T','N', nOcc(iSymJ), nOcc(iSymI),
     &                 nVir(iSymB), 2.0d0,
     &                 Wrk(kU+iOff),nVir(iSymB),
     &                 Wrk(kRic+iOff),nVir(iSymB),1.0d0,
     &                 Wrk(kPij(iSymJ)),nOcc(iSymJ))
 222              Continue
                  iOff1 = iOff1 + nOcc(iSymJ)*nVir(iSymB)
               End Do
            End Do

*           Calculate the "Exchange"-contribution to PiK (Froz-occ)
*           -------------------------------------------------------

            Do jVec1 = 1, NumLVecJ
                  iOffV1 = 0
                  iOffL1 = 0
                  Do iSymK = 1, nSym
                     iSymI = iSymK
                     iSymA = MulD2h(iSymK,jSym)
                     If(nOcc(iSymI)*nFro(iSymK)*nVir(iSymA) .eq. 0)
     &                      Go To 231
                     iOffV = nMoMo(jSym,iVecOV)*(jVec1-1) + iOffV1
                     iOffL = nMoMo(jSym,iVecFV)*(jVec1-1) + iOffL1
                     Call dGemm_('T','N', nOcc(iSymI), nFro(iSymK),
     &                          nVir(iSymA), 2.0d0,
     &                          Wrk(kV+iOffV),nVir(iSymA),
     &                          Wrk(kLKa+iOffL),nVir(iSymA),1.0d0,
     &                          Wrk(kPiK(iSymK)),nOcc(iSymI))
 231                 Continue
                     iOffV1 = iOffV1 + nOcc(iSymI)*nVir(iSymA)
                     iOffL1 = iOffL1 + nFro(iSymK)*nVir(iSymA)
                   End Do
               End Do

*           Calculate the "Exchange"-contribution to WiK (Froz-occ)
*           -------------------------------------------------------

            Do jVec1 = 1, NumLVecJ
                  iOffV1 = 0
                  iOffL1 = 0
                  Do iSymK = 1, nSym
                     iSymI = iSymK
                     iSymA = MulD2h(iSymK,jSym)
                     If(nOcc(iSymI)*nFro(iSymK)*nVir(iSymA) .eq. 0)
     &                      Go To 273
                     iOffV = nMoMo(jSym,iVecOV)*(jVec1-1) + iOffV1
                     iOffL = nMoMo(jSym,iVecFV)*(jVec1-1) + iOffL1
                     Call dGemm_('T','N', nOcc(iSymI), nFro(iSymK),
     &                          nVir(iSymA), -2.0d0,
     &                          Wrk(kV+iOffV),nVir(iSymA),
     &                          Wrk(kLKa+iOffL),nVir(iSymA),1.0d0,
     &                          Wrk(kWiK(iSymK)),nOcc(iSymI))
 273                 Continue
                     iOffV1 = iOffV1 + nOcc(iSymI)*nVir(iSymA)
                     iOffL1 = iOffL1 + nFro(iSymK)*nVir(iSymA)
                   End Do
               End Do


*           Calculate the "Exchange"-contribution to Lagr(IV)
*           -------------------------------------------------

            Do jVec1 = 1, NumLVecJ
               iOffV1 = 0
               iOffL1 = 0
               Do iSymI = 1, nSym
                  iSymA = iSymI
                  iSymB = MulD2h(iSymI,jSym)
                  If(nOcc(iSymI)*nVir(iSymB)*nVir(iSymA) .eq. 0)
     &                 Go To 223
                  iOffV = nMoMo(jSym,iVecOV)*(jVec1-1) + iOffV1
                  iOffL = nMoMo(jSym,iVecVV)*(jVec1-1) + iOffL1
                  Call dGemm_('T','N', nVir(iSymA), nOcc(iSymI),
     &                 nVir(iSymB), -2.0d0,
     &                 Wrk(kLab+iOffL),nVir(iSymB),
     &                 Wrk(kV+iOffV),nVir(iSymB),1.0d0,
     &                 Wrk(kLagr(iSymA)),nVir(iSymA))
 223              Continue
                  iOffL1 = iOffL1 + nVir(iSymA)*nVir(iSymB)
                  iOffV1 = iOffV1 + nOcc(iSymI)*nVir(iSymB)
               End Do
            End Do

*           Calculate the "Exchange"-contribution to Lagr(III)
*           --------------------------------------------------

            iOffV1 = 0
            iOffL1 = 0
            Do iSymJ = 1, nSym
               iSymA = MulD2h(jSym,iSymJ)
               iSymI = iSymA
               If(nOcc(iSymI)*nVir(iSymA)*nOcc(iSymJ).eq.0) Go To 233
               Do iJ = 1, nOcc(iSymJ)
                  iOffV = (iJ-1)*nVir(iSymA) + iOffV1
                  iOffL = (iJ-1)*nOcc(iSymI) + iOffL1
                  Call dGemm_('N','T',nVir(iSymA),nOcc(iSymI),
     &                       NumLVecJ, 2.0d0,
     &                       Wrk(kV+iOffV),nMoMo(jSym,iVecOV),
     &                       Wrk(kLji+iOffL),nMoMo(jSym,iVecOO),1.0d0,
     &                       Wrk(kLagr(iSymA)),nVir(iSymA))
               End Do
 233           Continue
               iOffV1 = iOffV1 + nOcc(iSymJ)*nVir(iSymA)
               iOffL1 = iOffL1 + nOcc(iSymJ)*nOcc(iSymI)
            End Do

*           Calculate the "Exchange"-contribution to Wai
*           --------------------------------------------------

            iOffV1 = 0
            iOffL1 = 0
            Do iSymJ = 1, nSym
               iSymA = MulD2h(jSym,iSymJ)
               iSymI = iSymA
               If(nOcc(iSymI)*nVir(iSymA)*nOcc(iSymJ)*NumLVecJ.eq.0)
     &                 Go To 2331
               Do iJ = 1, nOcc(iSymJ)
                  iOffV = (iJ-1)*nVir(iSymA) + iOffV1
                  iOffL = (iJ-1)*nOcc(iSymI) + iOffL1
                  Call dGemm_('N','T',nVir(iSymA),nOcc(iSymI),
     &                       NumLVecJ, -4.0d0,
     &                       Wrk(kV+iOffV),nMoMo(jSym,iVecOV),
     &                       Wrk(kLji+iOffL),nMoMo(jSym,iVecOO),1.0d0,
     &                       Wrk(kWai(iSymA)),nVir(iSymA))
               End Do
 2331          Continue
               iOffV1 = iOffV1 + nOcc(iSymJ)*nVir(iSymA)
               iOffL1 = iOffL1 + nOcc(iSymJ)*nOcc(iSymI)
            End Do

*           Calculate the "Exchange"-contribution to Wab
*           ---------------------------------------------
            iOff1 = 0
            Do iSymJ = 1, nSym
               iSymA = MulD2h(jSym,iSymJ)
               If(nOcc(iSymJ)*nVir(iSymA)*NumLVecJ.eq.0) Go To 251
               Do iI = 1, nOcc(iSymJ)
                  iOff = (iI-1)*nVir(iSymA) + iOff1
                  Call dGemm_('N','T',nVir(iSymA),nVir(iSymA),
     &                       NumLVecJ, -2.0d0,
     &                       Wrk(kLic+iOff),nMoMo(jSym,iVecOV),
     &                       Wrk(kV+iOff),nMoMo(jSym,iVecOV),1.0d0,
     &                       Wrk(kWab(iSymA)),nVir(iSymA))
               End Do
 251           Continue
               iOff1 = iOff1 + nOcc(iSymJ)*nVir(iSymA)
            End Do

*           Calculate the "Exchange"-contribution to Wij
*           --------------------------------------------

            Do jVec1 = 1, NumLVecJ
               iOff1 = 0
               Do iSymJ = 1, nSym
                  iSymI = iSymJ
                  iSymB = MulD2h(iSymJ,jSym)
                  If(nOcc(iSymJ)*nVir(iSymB) .eq. 0) Go To 252
                  iOff = nMoMo(jSym,iVecOV)*(jVec1-1) + iOff1
                  Call dGemm_('T','N', nOcc(iSymJ), nOcc(iSymI),
     &                 nVir(iSymB), -2.0d0,
     &                 Wrk(kV+iOff),nVir(iSymB),
     &                 Wrk(kLic+iOff),nVir(iSymB),1.0d0,
     &                 Wrk(kWij(iSymJ)),nOcc(iSymJ))
 252              Continue
                  iOff1 = iOff1 + nOcc(iSymJ)*nVir(iSymB)
               End Do
            End Do

*           Calculate "Exchange"-contribution to the Frozen Lagr(III)
*           --------------------------------------------------------
            iOffV1 = 0
            iOffL1 = 0
            Do iSymJ = 1, nSym
               iSymA = MulD2h(jSym,iSymJ)
               iSymI = iSymA
               If(nFro(iSymI)*nVir(iSymA)*nOcc(iSymJ).eq.0) Go To 2341
               Do iJ = 1, nOcc(iSymJ)
                  iOffV = (iJ-1)*nVir(iSymA) + iOffV1
                  iOffL = (iJ-1)*nFro(iSymI) + iOffL1
                  Call dGemm_('N','T',nVir(iSymA),nFro(iSymI),
     &                       NumLVecJ, 2.0d0,
     &                       Wrk(kV+iOffV),nMoMo(jSym,iVecOV),
     &                       Wrk(kLiK+iOffL),nMoMo(jSym,iVecOF),1.0d0,
     &                       Wrk(kFLagr(iSymA)),nVir(iSymA))
               End Do
 2341          Continue
               iOffV1 = iOffV1 + nOcc(iSymJ)*nVir(iSymA)
               iOffL1 = iOffL1 + nOcc(iSymJ)*nFro(iSymI)
            End Do

*           Calculate "Exchange"-contribution to the Wak(Vir-Froz)
*           --------------------------------------------------------
            iOffV1 = 0
            iOffL1 = 0
            Do iSymJ = 1, nSym
               iSymA = MulD2h(jSym,iSymJ)
               iSymI = iSymA
               If(nFro(iSymI)*nVir(iSymA)*nOcc(iSymJ).eq.0) Go To 234
               Do iJ = 1, nOcc(iSymJ)
                  iOffV = (iJ-1)*nVir(iSymA) + iOffV1
                  iOffL = (iJ-1)*nFro(iSymI) + iOffL1
                  Call dGemm_('N','T',nVir(iSymA),nFro(iSymI),
     &                       NumLVecJ, -4.0d0,
     &                       Wrk(kV+iOffV),nMoMo(jSym,iVecOV),
     &                       Wrk(kLiK+iOffL),nMoMo(jSym,iVecOF),1.0d0,
     &                       Wrk(kWaK(iSymA)),nVir(iSymA))
               End Do
 234           Continue
               iOffV1 = iOffV1 + nOcc(iSymJ)*nVir(iSymA)
               iOffL1 = iOffL1 + nOcc(iSymJ)*nFro(iSymI)
            End Do


         End Do

         Call DaClos(LuUVec)
         Call DaClos(LuVVec)
         Call ChoMP2_OpenF(iClos,iTypR,jSym)
         Call ChoMP2_OpenF(iClos,iTypL,jSym)
 10      Continue
      End Do                    ! jSym

      Call DaClos(LuRInv(1))
      Call DaClos(LuRInv(2))


*     Scale the Frozen-Occupied part of the density
*     ---------------------------------------------
      Do iSym = 1, nSym
         Do iI = 1, nOcc(iSym)
            E_i = EOcc(iOcc(iSym) + iI)
            Do iK = 1, nFro(iSym)
               E_K = EFro(iFro(iSym) + iK)
               i = iI-1 + (iK-1)*nOcc(iSym)
               Wrk(kPik(iSym)+i) = Wrk(kPik(iSym)+i)/(E_i - E_k)
            End Do
         End Do
      End Do



      Do iSym = 1, nSym
         If(NumCho(iSym) .gt. 0) Then

*
            Call ChoMP2_OpenF(1,1,iSym)

*
*           I doubt the intermediate files here need to be
*           Multifiles, probably they will never be used again after the
*           end of the iSym-loop (or even IB-block loop).
            iSeed = 7
            LuUVec = IsFreeUnit(iSeed)
            Write(Fname,'(A4,I1)') 'TMPV',1
            Call DaName_MF_WA(LuUVec,Fname)
*
            iSeed = 7
            LuVVec = IsFreeUnit(iSeed)
            Write(Fname,'(A4,I1)') 'TMPV',2
            Call DaName_MF_WA(LuVVec,Fname)
*
            iSeed = 7
            LuWVec = IsFreeUnit(iSeed)
            Write(Fname,'(A4,I1)') 'TMPV',3
            Call DaName_MF_WA(LuWVec,Fname)
*
            lScr = lWrk - kEndFLagr
*           Construct Lagr(i)
*           -----------------
            Call ChoMP2g_ConstrAP(irc,Wrk(kEndFLagr),lScr,
     &                            'oovo',iSym,nVec,Wrk(kLagr(1)),
     &                            lLagr,Wrk(kPij(1)),lPij,-0.5d0)

            Call ChoMP2g_ConstrAP(irc,Wrk(kEndFLagr),lScr,
     &                            'fovo',iSym,nVec,Wrk(kLagr(1)),
     &                            lLagr,Wrk(kPiK(1)),lPiK,-1.0d0)
*           Construct FLagr(i)
*           ------------------
            Call ChoMP2g_ConstrAP(irc,Wrk(kEndFLagr),lScr,
     &                            'oovf',iSym,nVec,Wrk(kFLagr(1)),
     &                            lFLagr,Wrk(kPij(1)),lPij,-0.5d0)

            Call ChoMP2g_ConstrAP(irc,Wrk(kEndFLagr),lScr,
     &                            'fovf',iSym,nVec,Wrk(kFLagr(1)),
     &                            lFLagr,Wrk(kPiK(1)),lPiK,-1.0d0)
*           Construct Lagr(ii)
*           ------------------
            Call ChoMP2g_ConstrAP(irc,Wrk(kEndFLagr),lScr,
     &                            'vvvo',iSym,nVec,Wrk(kLagr(1)),
     &                            lLagr,Wrk(kPab(1)),lPab,-0.5d0)
*           Construct FLagr(ii)
*           -------------------
            Call ChoMP2g_ConstrAP(irc,Wrk(kEndFLagr),lScr,
     &                            'vvvf',iSym,nVec,Wrk(kFLagr(1)),
     &                            lFLagr,Wrk(kPab(1)),lPab,-0.5d0)

            Call ChoMP2_OpenF(iClos,1,iSym)
            Call DaClos(LuUVec)
            Call DaClos(LuVVec)
            Call DaClos(LuWVec)

         End If
      End Do
*

*
#ifdef _DEBUGPRINT_
      Write(6,*) 'Pab'
      Do i = 1,lPab
         Write(6,*) Wrk(kPab(1)+i-1)
      End Do
      Write(6,*) 'lWab'
      Do i = 1, lWab
         Write(6,*) Wrk(kWab(1)+i-1)
      End Do
      Write(6,*) 'Pij'
      Do i = 1, lPij
         Write(6,*) Wrk(kPij(1)+i-1)
      End Do
      Write(6,*) 'Wij'
      Do i = 1, lWij
         Write(6,*) Wrk(kWij(1)+i-1)
      End Do
      Write(6,*) 'WiK'
      Do i = 1, lWiK
         Write(6,*) Wrk(kWiK(1)+i-1)
      End Do
      Write(6,*) 'WaK'
      Do i = 1, lWaK
         Write(6,*) Wrk(kWaK(1)+i-1)
      End Do
      Write(6,*) 'Wai'
      Do i = 1, lWai
         Write(6,*) Wrk(kWai(1)+i-1)
      End Do
      Write(6,*) 'PiK'
      Do i = 1,lPiK
         Write(6,*) Wrk(kPiK(1)+i-1)
      End Do

      Write(6,*) 'Lagr'
      Do i = 1,lLagr
         Write(6,*) Wrk(kLagr(1)+i-1)
      End Do

      Write(6,*) 'FLagr'
      Do i = 1,lFLagr
         Write(6,*) Wrk(kFLagr(1)+i-1)
      End Do
#endif
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(EVir)
      End
