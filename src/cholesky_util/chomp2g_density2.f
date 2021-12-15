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

      SubRoutine ChoMP2g_density2(irc,EOcc,EVir,EFro,Wrk,lWrk)
*
*     Jonas Bostrom, Mars 2010
*
*     Purpose: Solve the CPHF-equations to obtain occ-vir
*              contributions to the 1-pdm.

#include "implicit.fh"
#include "chomp2g.fh"
#include "chomp2.fh"
#include "cholesky.fh"
#include "choorb.fh"
#include "WrkSpc.fh"
*
      Character Fname*5
      Real*8 Wrk(lWrk),EOcc(*), EVir(*), EFro(*)
      Integer kCGVec(9), kEndCGVec(9), kDiag(2), nOccAll(8)
      Logical Done

      Character*8  ThisNm
      Character*16 SecNam
      Parameter (SecNam = 'ChoMP2g_Density2', ThisNm = 'Density2')

*     Do not delete vectors after use.
*     ------------------------------------------------------
      iClos = 2
      iTypL = 1
      iVecFF = 1
      iVecOV = 6
      iVecVV = 9
      iVecOO = 5
      iVecFV = 3
      iVecFO = 2
      iVecOF = 4

*
*     Setup the conjugate gradient procedure
*     --------------------------------------
      Eps = 1D-8
      Done=.false.
      nIter = 100

*     Allocate memory for Frozen Diagonal of A
*     ---------------------------------
      lFDiag = nMoMo(1,iVecFV)
      kDiag(1) = kFLagr(1) + lFLagr
      kEndFDiag = kDiag(1) + lFDiag
      Call FZero(Wrk(kDiag(1)),lFDiag)

*     Allocate memory for Occupied Diagonal of A
*     ---------------------------------
      lDiag = nMoMo(1,iVecOV)
      kDiag(2) = kEndFDiag
      kEndDiag = kDiag(2) + lDiag
      Call FZero(Wrk(kDiag(2)),lDiag)
*
      iSym = 1
      nOccAll(iSym) = nFro(iSym)+nOcc(iSym)

*     Open Cholesky vector files.
*     ---------------------------
      Call ChoMP2_OpenF(1,1,iSym)

      maxvalue = 200
      nVec = min(NumCho(iSym),maxvalue)
      nBatL = (NumCho(iSym)-1)/nVec + 1
      If(nVec .lt. 1) Then
         Call ChoMP2_Quit(SecNam,'Insufficient memory','[1]')
      End If

*     Allocate memory for Lia-vector and LIa-vector
*     ---------------------------------------------

      lLfa = nMoMo(iSym,iVecFV)*nVec
      kLfa = kEndDiag
      kEndLfa = kLfa + lLfa

      lLia = nMoMo(iSym,iVecOV)*nVec
      kLia = kEndLfa
      kEndLia = kLia + lLia

      Do iBat = 1, nBatL
         If (iBat .eq. nBatL) Then
            NumVec = NumCho(iSym) - nVec*(nBatL-1)
         Else
            NumVec = nVec
         End If
         iVec = nVec*(iBat-1) + 1

*        Read Lfa-vectors
*        ----------------
         iOpt = 2
         lTot = nMoMo(iSym,iVecFV)*NumVec
         iAdr = nMoMo(iSym,iVecFV)*(iVec-1) + 1 +
     &        iAdrOff(iSym,iVecFV)
         Call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLfa),
     &                lTot,iAdr)

*        Read Lia-vectors
*        ----------------
         iOpt = 2
         lTot = nMoMo(iSym,iVecOV)*NumVec
         iAdr = nMoMo(iSym,iVecOV)*(iVec-1) + 1 +
     &        iAdrOff(iSym,iVecOV)
         Call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLia),
     &                lTot,iAdr)

*        Construct Diagonal of A
*        -----------------------
         Do iVec1 = 1, NumVec
            Do i = 1, nMoMo(iSym,iVecFV)
               iOffL = (iVec1-1)*nMoMo(iSym,iVecFV)
               Wrk(kDiag(1)+i-1) =  Wrk(kDiag(1)+i-1) +
     &                           3.0d0*Wrk(kLfa+i-1+iOffL)**2
            End Do
         End Do
         Do iVec1 = 1, NumVec
            Do i = 1, nMoMo(iSym,iVecOV)
               iOffL = (iVec1-1)*nMoMo(iSym,iVecOV)
               Wrk(kDiag(2)+i-1) = Wrk(kDiag(2)+i-1) +
     &                             3.0d0*Wrk(kLia+i-1+iOffL)**2
            End Do
         End Do
*
         index1 = 0
         Do iSym1 = 1, nSym
            nA = nVir(iSym1)
            nI = nFro(iSym1)
            Do iI = 1, nI
               Ei = EFro(iFro(iSym1) + iI)
               Do iA = 1,nA
                  index = iA-1 + (iI-1)*nA + index1
                  Ea = EVir(iVir(iSym1) + iA)
                  Wrk(kDiag(1)+index) = Wrk(kDiag(1)+index) + Ea - Ei
               End Do
            End Do
            nI = nOcc(iSym1)
            Do iI = 1, nI
               Ei = EOcc(iOcc(iSym1) + iI)
               Do iA = 1,nA
                  index = iA-1 + (iI-1)*nA + index1
                  Ea = EVir(iVir(iSym1) + iA)
                  Wrk(kDiag(2)+index) = Wrk(kDiag(2)+index) + Ea - Ei
               End Do
            End Do
            index1 = index1 + nOcc(iSym1)*nVir(iSym1)
         End Do
         Do i = 1, lFDiag
            Wrk(kDiag(1)+i-1) = 1/Wrk(kDiag(1)+i-1)
         End Do
         Do i = 1, lDiag
            Wrk(kDiag(2)+i-1) = 1/Wrk(kDiag(2)+i-1)
         End Do
      End Do

      Call ChoMP2_OpenF(iClos,1,iSym)


*     Allocate vectors needed for the PCG
*     -----------------------------------
      lCGFVec = 0
      Do iSym = 1, nSym
         lCGFVec = lCGFVec + (nFro(iSym))*nVir(iSym)
      End Do
      lCGOVec = 0
      Do iSym = 1, nSym
         lCGOVec = lCGOVec + (nOcc(iSym))*nVir(iSym)
      End Do
      lCGVec = lCGOVec + lCGFVec
*     Vector Legend (as they are named in Conj_Grad):
*                    Z = 1, Ztemp = 2
*                    R = 3, Rtemp = 4
*                    P = 5, Ptemp = 6
*                    X = 7, Xtemp = 8
*                    AP = 9

      kCGVec(1) = kEndDiag
      kEndCGVec(1) = kCGVec(1) + lCGVec
      Call FZero(Wrk(kCGVec(1)),lCGVec)
      nCGvec = 9

      Do i = 2, nCGVec
         kCGVec(i) = kEndCGVec(i-1)
         kEndCGVec(i) = kCGVec(i) + lCGVec
         Call FZero(Wrk(kCGVec(i)),lCGVec)
      End Do

*     Calculate inital values for the CG-vectors needing that
*     -------------------------------------------------------

      Do i=1, lCGFVec
         Wrk(kCGVec(1)+i-1) = Wrk(kFLagr(1)+i-1)*Wrk(kDiag(1)+i-1)
         Wrk(kCGVec(3)+i-1) = Wrk(kFLagr(1)+i-1)
         Wrk(kCGVec(5)+i-1) = Wrk(kCGVec(1)+i-1)
      End Do
      iOff = lCGFVec
      Do i = 1, lCGOVec
         Wrk(kCGVec(1)+i-1+iOff) = Wrk(kLagr(1)+i-1)*Wrk(kDiag(2)+i-1)
         Wrk(kCGVec(3)+i-1+iOff) = Wrk(kLagr(1)+i-1)
         Wrk(kCGVec(5)+i-1+iOff) = Wrk(kCGVec(1)+i-1+iOff)
      End Do

      Do iIter = 1, nIter
         Call FZero(Wrk(kCGVec(9)),lCGVec)

*     Calculate A*p
*     -------------
         Do iSym = 1, nSym

*           Open Cholesky vector files.
*           ---------------------------
            Call ChoMP2_OpenF(1,1,iSym)

            iSeed = 7
            LuVVec = IsFreeUnit(iSeed)
            Write(Fname,'(A4,I1)') 'TMPV',2
            Call DaName_MF_WA(LuVVec,Fname)
*
            iSeed = 7
            LuWVec = IsFreeUnit(iSeed)
            Write(Fname,'(A4,I1)') 'TMPV',3
            Call DaName_MF_WA(LuWVec,Fname)

            nVec = min(NumCho(iSym),maxvalue)
            If(nVec .lt. 1) Then
               Call ChoMP2_Quit(SecNam,'Insufficient memory','[1]')
            End If

            lScr = lWrk - kEndCGVec(9)
            iOff = lCGFVec
*
*           Construct The Frozen part of A*P
*           --------------------------------
            Call ChoMP2g_Constrap(irc,Wrk(kEndCGVec(9)),lScr,
     &                            'fvvf',iSym,nVec,Wrk(kCGVec(9)),
     &                            lCGFVec,Wrk(kCGVec(5)),lCGFVec,
     &                            1.0d0)
            Call ChoMP2g_Constrap(irc,Wrk(kEndCGVec(9)),lScr,
     &                            'ovvf',iSym,nVec,Wrk(kCGVec(9)),
     &                            lCGFVec,Wrk(kCGVec(5)+iOff),lCGOVec,
     &                            1.0d0)
*           Construct The Occupied part of A*P
*           ----------------------------------
            Call ChoMP2g_Constrap(irc,Wrk(kEndCGVec(9)),lScr,
     &                            'fvvo',iSym,nVec,Wrk(kCGVec(9)+iOff),
     &                            lCGOVec,Wrk(kCGVec(5)),lCGFVec,
     &                            1.0d0)
            Call ChoMP2g_Constrap(irc,Wrk(kEndCGVec(9)),lScr,
     &                            'ovvo',iSym,nVec,Wrk(kCGVec(9)+iOff),
     &                            lCGOVec,Wrk(kCGVec(5)+iOff),lCGOVec,
     &                            1.0d0)

            Call ChoMP2_OpenF(iClos,1,iSym)
            Call DaClos(LuVVec)
            Call DaClos(LuWVec)

         End Do
*
         index1 = 0
         Do iSym1 = 1, nSym
            nA = nVir(iSym1)
            nI = nFro(iSym1)
            Do iI = 1, nI
               Ei = EFro(iFro(iSym1) + iI)
               Do iA = 1,nA
                  index = iA-1 + (iI-1)*nA + index1
                  Ea = EVir(iVir(iSym1) + iA)
                  Wrk(kCGVec(9)+ index) =Wrk(kCGVec(9)+index)+
     &                                   (Ea - Ei)*Wrk(kCGVec(5)+index)
               End Do
            End Do
            nI = nOcc(iSym1)
            Do iI = 1, nI
               Ei = EOcc(iOcc(iSym1) + iI)
               Do iA = 1,nA
                  index = iA-1 + (iI-1)*nA +index1
                  Ea = EVir(iVir(iSym1) + iA)
                  Wrk(kCGVec(9)+ index+iOff) =Wrk(kCGVec(9)+index+iOff)+
     &                              (Ea - Ei)*Wrk(kCGVec(5)+index+iOff)
               End Do
            End Do
            index1 = index1 + (nFro(iSym1)+nOcc(iSym1))*nVir(iSym1)
         End Do
*
         Call Conj_Grad(Done,lCGVec,Wrk(kDiag(1)),Wrk(kCGVec(7)),
     &                  Wrk(kCGVec(8)),Wrk(kCGVec(3)),Wrk(kCGVec(4)),
     &                  Wrk(kCGVec(5)),Wrk(kCGVec(6)),Wrk(kCGVec(1)),
     &                  Wrk(kCGVec(2)),Wrk(kCGVec(9)),Eps,Res)
         If(Done) Go To 100
      End Do
      Write(6,*) '***************WARNING************************'
      Write(6,*) ''
      write(6,*) 'Too many iterations, this is what you get after:'
      Write(6,*)  nIter, ' Iterations'
      write(6,*) 'The residual is', res, 'and not', Eps
      Write(6,*) '**********************************************'

 100  Continue

*     Construct MP2 Density contribution from parts of the matrix
*     -----------------------------------------------------------
      iSymOffOV = 0
      Do iSym = 1, nSym
         Do i = 1, nOcc(iSym)
            iOrb = nFro(iSym) + i
            Do j = 1, nFro(iSym)
               jOrb = j
               Work(ipDensity(iSym) + jOrb-1 + (iOrb-1)*nOrb(iSym)) =
     &             Wrk(kPiK(iSym) + i-1 + nOcc(iSym)*(j-1))
               Work(ipDensity(iSym) + iOrb-1 + (jOrb-1)*nOrb(iSym)) =
     &             Wrk(kPiK(iSym) + i-1 + nOcc(iSym)*(j-1))
            End Do
         End Do

         Do i = 1, nOcc(iSym)
            iOrb = nFro(iSym) + i
            Do j =  1, nOcc(iSym)
               jOrb = nFro(iSym) + j
               Work(ipDensity(iSym) + jOrb-1 + (iOrb-1)*nOrb(iSym)) =
     &               Wrk(kPij(iSym) + j-1 + nOcc(iSym)*(i-1))
            End Do
         End  Do

         Do iI = 1, nFro(iSym)
            Do iA = 1, nVir(iSym)
               Wrk(kPaK(iSym)+ iA-1 + nVir(iSym)*(iI-1)) =
     &           Wrk(kCGVec(7) + iA-1 + nVir(iSym)*(iI-1) + iSymOffOV)
            End Do
         End Do
         Do iI = 1, nOcc(iSym)
            iOrbI = nFro(iSym) + iI
            Do iA = 1, nVir(iSym)
               Wrk(kPai(iSym)+ iA-1 + nVir(iSym)*(iI-1)) =
     &          Wrk(kCGVec(7) + iA-1 + nVir(iSym)*(iOrbI-1) + iSymOffOV)
            End Do
         End Do

         Do i = 1, nFro(iSym) + nOcc(iSym)
            iOrb = i
            Do j = 1, nVir(iSym)
               jOrb = nFro(iSym) + nOcc(iSym) +  j
               Work(ipDensity(iSym) + jOrb-1 + (iOrb-1)*nOrb(iSym)) =
     &               Wrk(kCGVec(7) + j-1 +
     &               (nVir(iSym))*(i-1) + iSymOffOV)
               Work(ipDensity(iSym) + iOrb-1 + (jOrb-1)*nOrb(iSym)) =
     &             Work(ipDensity(iSym) + jOrb-1 + (iOrb-1)*nOrb(iSym))

            End Do
         End Do

         Do i = 1, nVir(iSym)
            iOrb = nFro(iSym) + nOcc(iSym) + i
            Do j =  1, nVir(iSym)
               jOrb = nFro(iSym) + nOcc(iSym) + j
               Work(ipDensity(iSym) + jOrb-1 + (iOrb-1)*nOrb(iSym)) =
     &               Wrk(kPab(iSym) + j-1 + nVir(iSym)*(i-1))
            End Do
         End  Do
         iSymOffOV = iSymOffOV + (nOcc(iSym)+nFro(iSym))*nVir(iSym)
      End Do

*     Add type (II) terms to W-matrix
*     -------------------------------

      Do iSym = 1, nsym
         Do iI = 1, nOcc(iSym)
            Ei = EOcc(iOcc(iSym)+iI)
            Do iJ = 1, nFro(iSym)
               Ej = EFro(iFro(iSym)+iJ)
               Wrk(kWiK(iSym) + iI-1 + nOcc(iSym)*(iJ-1)) =
     &             Wrk(kWiK(iSym) + iI-1 + nOcc(iSym)*(iJ-1)) -1.0*
     &             Wrk(kPiK(iSym) + iI-1 + nOcc(iSym)*(iJ-1))*(Ei+Ej)
            End Do
         End Do
         Do iI = 1, nOcc(iSym)
            Ei = EOcc(iOcc(iSym)+iI)
            Do iJ = 1, nOcc(iSym)
               Ej = EOcc(iOcc(iSym)+iJ)
               Wrk(kWij(iSym) + iJ-1 + nOcc(iSym)*(iI-1)) =
     &             Wrk(kWij(iSym) + iJ-1 + nOcc(iSym)*(iI-1)) - 0.5d0*
     &             Wrk(kPij(iSym) + iJ-1 + nOcc(iSym)*(iI-1))*(Ei+Ej)
            End Do
         End Do
         Do iA = 1, nVir(iSym)
            Ea = EVir(iVir(iSym)+iA)
            Do iB = 1, nVir(iSym)
               Eb = EVir(iVir(iSym)+iB)
               Wrk(kWab(iSym) + iB-1 + nVir(iSym)*(iA-1)) =
     &              Wrk(kWab(iSym) + iB-1 + nVir(iSym)*(iA-1)) - 0.5d0*
     &              Wrk(kPab(iSym) + iB-1 + nVir(iSym)*(iA-1))*(Ea+Eb)
            End Do
         End Do
         Do iI = 1, nFro(iSym)
            Ei = Efro(iFro(iSym)+iI)
            iOrbI = iI
            Do iA = 1, nVir(iSym)
               iOrbA = nFro(iSym) + nOcc(iSym) + iA
               Wrk(kWak(iSym) + iA-1 + nVir(iSym)*(iI-1)) =
     &           Wrk(kWak(iSym) + iA-1 + nVir(iSym)*(iI-1)) - 2.0d0*
     &           Wrk(kPaK(iSym) + iA-1 + nVir(iSym)*(iI-1))*Ei
            End Do
         End Do
         Do iI = 1, nOcc(iSym)
            Ei = EOcc(iOcc(iSym)+iI)
            iOrbI = nFro(iSYm) + iI
            Do iA = 1, nVir(iSym)
               iOrbA = nFro(iSym) + nOcc(iSym) + iA
               Wrk(kWai(iSym) + iA-1 + nVir(iSym)*(iI-1)) =
     &           Wrk(kWai(iSym) + iA-1 + nVir(iSym)*(iI-1)) - 2.0d0*
     &           Wrk(kPai(iSym) + iA-1 + nVir(iSym)*(iI-1))*Ei
            End Do
         End Do

      End Do

      iSym = 1

*     Add type (III) terms to W-matrix
*     -------------------------------

*     Open Cholesky vector files.
*     ---------------------------
      Call ChoMP2_OpenF(1,1,iSym)


      nVec = min(NumCho(iSym),maxvalue)
      If(nVec .lt. 1) Then
         Call ChoMP2_Quit(SecNam,'Insufficient memory','[1]')
      End If
      nBatL = (NumCho(iSym)-1)/nVec + 1

*     Allocate memory for U-vector
*     ----------------------------

      lU = nVec
      kU = kEndDiag
      kEndU = kU + lU


*     Allocate memory for Lia-vector and LIa-vector
*     ---------------------------------------------

      lLfa = nMoMo(iSym,iVecFV)*nVec
      kLfa = kEndU
      kEndLfa = kLfa + lLfa

      lLia = nMoMo(iSym,iVecOV)*nVec
      kLia = kEndLfa
      kEndLia = kLia + lLia

      lLij = nMoMo(iSym,iVecOO)*nVec
      kLij = kEndLia
      kEndLij = kLij + lLij

      lLiK = nMoMo(iSym,iVecOF)*nVec
      kLiK = kEndLij
      kEndLiK = kLiK + lLiK

      lLKi = nMoMo(iSym,iVecFO)*nVec
      kLKi = kEndLiK
      kEndLKi = kLKi + lLKi

      lLab = nMoMo(iSym,iVecVV)*nVec
      kLab = kEndLKi
      kEndLab = kLab + lLab

      lLJK = nMoMo(iSym,iVecFF)*nVec
      kLJK = kEndLab
      kEndLJK = kLJK + lLJK

      Do iBat = 1, nBatL
         If (iBat .eq. nBatL) Then
            NumVec = NumCho(iSym) - nVec*(nBatL-1)
         Else
            NumVec = nVec
         End If
         iVec = nVec*(iBat-1) + 1
         Call FZero(Wrk(kU),lU)

*     Read Lij^J-vectors from disk
*     ----------------------------

         If(NumVec .gt. 0) Then
            iOpt = 2
            lTot = nMoMo(iSym,iVecOO)*NumVec
            iAdr = 1 + nMoMo(iSym,iVecOO)*(iVec-1) +
     &           iAdrOff(iSym,iVecOO)
            Call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLij),lTot,
     &           iAdr)
         End If

*     Read LiK^J-vectors from disk
*     ----------------------------

         If(NumVec .gt. 0) Then
            iOpt = 2
            lTot = nMoMo(iSym,iVecOF)*NumVec
            iAdr = 1 + nMoMo(iSym,iVecOF)*(iVec-1) +
     &           iAdrOff(iSym,iVecOF)
            Call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLiK),lTot,
     &           iAdr)
         End If

*     Read LKi^J-vectors from disk
*     ----------------------------

         If(NumVec .gt. 0) Then
            iOpt = 2
            lTot = nMoMo(iSym,iVecFO)*NumVec
            iAdr = 1 + nMoMo(iSym,iVecFO)*(iVec-1) +
     &           iAdrOff(iSym,iVecFO)
            Call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLKi),lTot,
     &           iAdr)
         End If

*     Read LIK^J-vectors from disk
*     ----------------------------

         If(NumVec .gt. 0) Then
            iOpt = 2
            lTot = nMoMo(iSym,iVecFF)*NumVec
            iAdr = 1 + nMoMo(iSym,iVecFF)*(iVec-1) +
     &           iAdrOff(iSym,iVecFF)
            Call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLJK),lTot,
     &           iAdr)
         End If

*     Read Lia^J-vectors from disk
*     ----------------------------

         If(NumVec .gt. 0) Then
            iOpt = 2
            lTot = nMoMo(iSym,iVecOV)*NumVec
            iAdr = 1 + nMoMo(iSym,iVecOV)*(iVec-1) +
     &           iAdrOff(iSym,iVecOV)
            Call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLia),lTot,
     &           iAdr)
         End If

*     Read Lfa^J-vectors from disk
*     ----------------------------

         If(NumVec .gt. 0) Then
            iOpt = 2
            lTot = nMoMo(iSym,iVecFV)*NumVec
            iAdr = 1 + nMoMo(iSym,iVecFV)*(iVec-1) +
     &           iAdrOff(iSym,iVecFV)
            Call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLfa),lTot,
     &           iAdr)
         End If

*     Read Lab^J-vectors from disk
*     ----------------------------

         If(NumVec .gt. 0) Then
            iOpt = 2
            lTot = nMoMo(iSym,iVecVV)*NumVec
            iAdr = 1 + nMoMo(iSym,iVecVV)*(iVec-1) +
     &           iAdrOff(iSym,iVecVV)
            Call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLab),lTot,
     &           iAdr)
         End If

*     Construct U^J intermediate vectors
*     ----------------------------------
         If(NumVec*nMoMo(iSym,iVecOO).eq.0) Go To 101
         Call dGemm_('N','N',1,NumVec,nMoMo(iSym,iVecOO),
     &              1.0d0,Wrk(kPij(iSym)), 1,
     &              Wrk(kLij),nMoMo(iSym,iVecOO),1.0d0,
     &              Wrk(kU),1)
 101     Continue

         If(NumVec*nMoMo(iSym,iVecFO).eq.0) Go To 102
         Call dGemm_('N','N',1,NumVec,nMoMo(iSym,iVecFO),
     &              2.0d0,Wrk(kPiK(iSym)), 1,
     &              Wrk(kLKi),nMoMo(iSym,iVecFO),1.0d0,
     &              Wrk(kU),1)
 102     Continue

         If(NumVec*nMoMo(iSym,iVecOV).eq.0) Go To 103
         Call dGemm_('N','N',1,NumVec,nMoMo(iSym,iVecOV),
     &              2.0d0,Wrk(kPai(iSym)), 1,
     &              Wrk(kLia),nMoMo(iSym,iVecOV),1.0d0,
     &              Wrk(kU),1)
 103     Continue

         If(NumVec*nMoMo(iSym,iVecFV).eq.0) Go To 104
         Call dGemm_('N','N',1,NumVec,nMoMo(iSym,iVecFV),
     &              2.0d0,Wrk(kPaK(iSym)), 1,
     &              Wrk(kLfa),nMoMo(iSym,iVecFV),1.0d0,
     &              Wrk(kU),1)
 104     Continue

         If(NumVec*nMoMo(iSym,iVecVV).eq.0) Go To 105
         Call dGemm_('N','N',1,NumVec,nMoMo(iSym,iVecVV),
     &              1.0d0,Wrk(kPab(iSym)), 1,
     &              Wrk(kLab),nMoMo(iSym,iVecVV),1.0d0,
     &              Wrk(kU),1)
 105     Continue

*     Construct contribution to Wij
*     -----------------------------

         If(nMoMo(iSym,iVecOO).eq.0) Go To 111
         Call dGemm_('N','N',nMoMo(iSym,iVecOO),1,NumVec,
     &             -2.0d0,Wrk(kLij),nMoMo(iSym,iVecOO),
     &             Wrk(kU),NumVec,1.0d0,
     &             Wrk(kWij(iSym)),nMoMo(iSym,iVecOO))
 111     Continue

*     Construct Contribution to WiK
*     -----------------------------

         If(nMoMo(iSym,iVecOF).eq.0) Go To 112
         Call dGemm_('N','N',nMoMo(iSym,iVecOF),1,NumVec,
     &             -4.0d0,Wrk(kLKi),nMoMo(iSym,iVecOF),
     &             Wrk(kU),NumVec,1.0d0,
     &             Wrk(kWiK(iSym)),nMoMo(iSym,iVecOF))
 112     Continue

*     Construct Contribution to WJK
*     -----------------------------

         If(nMoMo(iSym,iVecFF).eq.0) Go To 113
         Call dGemm_('N','N',nMoMo(iSym,iVecFF),1,NumVec,
     &             -2.0d0,Wrk(kLJK),nMoMo(iSym,iVecFF),
     &             Wrk(kU),NumVec,1.0d0,
     &             Wrk(kWJK(iSym)),nMoMo(iSym,iVecFF))
 113     Continue

      End Do


      iVecFF = 1
      iVecOF = 2
      iVecVF = 3
      iVecFO = 4
      iVecOO = 5
      iVecVO = 6
      iVecFV = 7
      iVecOV = 8
      iVecVV = 9

*     Allocate memory for Lia-vector and LIa-vector
*     ---------------------------------------------

      lLJK = nMoMo(iSym,iVecFF)*nVec
      kLJK = kEndDiag
      kEndLJK = kLJK + lLJK

      lLKi = nMoMo(iSym,iVecOF)*nVec
      kLKi = kEndLJK
      kEndLKi = kLKi + lLKi

      lLKa = nMoMo(iSym,iVecVF)*nVec
      kLKa = kEndLKi
      kEndLKa = kLKa + lLKa

      lLiK = nMoMo(iSym,iVecFO)*nVec
      kLiK = kEndLKa
      kEndLiK = kLiK + lLiK

      lLij = nMoMo(iSym,iVecOO)*nVec
      kLij = kEndLiK
      kEndLij = kLij + lLij

      lLia = nMoMo(iSym,iVecVO)*nVec
      kLia = kEndLij
      kEndLia = kLia + lLia

      lLip = nOccAll(iSym)*nOrb(iSym)*nVec
      kLip = kEndLia
      kEndLip = kLip + lLip

      lVip = nOccAll(iSym)*nOrb(iSym)*nVec
      kVip = kEndLip
      kEndVip = kVip + lVip

      lWij2 = nOccAll(iSym)*nOccAll(iSym)
      kWij2 = kEndVip
      kEndWij2 = kWij2 + lWij2
      Do i = 1, lWij2
         Wrk(kWij2+i-1) = 0.0d0
      End Do

      Do iBat = 1, nBatL
         If (iBat .eq. nBatL) Then
            NumVec = NumCho(iSym) - nVec*(nBatL-1)
         Else
            NumVec = nVec
         End If
         iVec = nVec*(iBat-1) + 1

*     Read Lpq^J-vectors from disk
*     ----------------------------

         If(NumVec .gt. 0) Then
            iOpt = 2
            lTot = nMoMo(iSym,iVecFF)*NumVec
            iAdr = 1 + nMoMo(iSym,iVecFF)*(iVec-1) +
     &           iAdrOff(iSym,iVecFF)
            Call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLJK),lTot,
     &           iAdr)

            lTot = nMoMo(iSym,iVecOF)*NumVec
            iAdr = 1 + nMoMo(iSym,iVecOF)*(iVec-1) +
     &           iAdrOff(iSym,iVecOF)
            Call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLKi),lTot,
     &           iAdr)

            lTot = nMoMo(iSym,iVecVF)*NumVec
            iAdr = 1 + nMoMo(iSym,iVecVF)*(iVec-1) +
     &           iAdrOff(iSym,iVecVF)
            Call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLKa),lTot,
     &           iAdr)

            lTot = nMoMo(iSym,iVecFO)*NumVec
            iAdr = 1 + nMoMo(iSym,iVecFO)*(iVec-1) +
     &           iAdrOff(iSym,iVecFO)
            Call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLiK),lTot,
     &           iAdr)

            lTot = nMoMo(iSym,iVecOO)*NumVec
            iAdr = 1 + nMoMo(iSym,iVecOO)*(iVec-1) +
     &           iAdrOff(iSym,iVecOO)
            Call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLij),lTot,
     &           iAdr)

            lTot = nMoMo(iSym,iVecVO)*NumVec
            iAdr = 1 + nMoMo(iSym,iVecVO)*(iVec-1) +
     &           iAdrOff(iSym,iVecVO)
            Call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLia),lTot,
     &           iAdr)
         End If

*     Put the Lpq-vectors together to one large vector
*     ------------------------------------------------

         iOff2 = 0
         Do iVec1 = 1, NumVec
            Do i = 1, nFro(iSym)
               iOff1 = nMoMo(iSym,iVecFF)*(iVec1-1) + (i-1)*nFro(iSym)
               Call dCopy_(nFro(iSym),Wrk(kLJK+iOff1),1,
     &                               Wrk(kLip+iOff2),1)
               iOff2 = iOff2+nFro(iSym)
               iOff1 = nMoMo(iSym,iVecOF)*(iVec1-1) + (i-1)*nOcc(iSym)
               Call dCopy_(nOcc(iSym),Wrk(kLKi+iOff1),1,
     &                               Wrk(kLip+iOff2),1)
               iOff2 = iOff2+nOcc(iSym)
               iOff1 = nMoMo(iSym,iVecVF)*(iVec1-1) + (i-1)*nVir(iSym)
               Call dCopy_(nVir(iSym),Wrk(kLKa+iOff1),1,
     &                               Wrk(kLip+iOff2),1)
               iOff2 = iOff2+nVir(iSym)
            End Do
            Do i = 1, nOcc(iSym)
               iOff1 = nMoMo(iSym,iVecFO)*(iVec1-1) + (i-1)*nFro(iSym)
               Call dCopy_(nFro(iSym),Wrk(kLiK+iOff1),1,
     &                               Wrk(kLip+iOff2),1)
               iOff2 = iOff2 + nFro(iSym)
               iOff1 = nMoMo(iSym,iVecOO)*(iVec1-1) + (i-1)*nOcc(iSym)
               Call dCopy_(nOcc(iSym),Wrk(kLij+iOff1),1,
     &                               Wrk(kLip+iOff2),1)
               iOff2 = iOff2 + nOcc(iSym)
               iOff1 = nMoMo(iSym,iVecVO)*(iVec1-1) + (i-1)*nVir(iSym)
               Call dCopy_(nVir(iSym),Wrk(kLia+iOff1),1,
     &                               Wrk(kLip+iOff2),1)
               iOff2 = iOff2 + nVir(iSym)
            End Do
         End Do

         Do iVec1 = 1, NumVec
            iOff = nOrb(iSym)*nOccAll(iSym)*(iVec1-1)
            Call dGemm_('N','N',nOrb(iSym),nOccAll(iSym),nOrb(iSym),
     &                  1.0d0,Work(ipDensity(iSym)),nOrb(iSym),
     &                  Wrk(kLip+iOff),nOrb(iSym),0.0d0,
     &                  Wrk(kVip+iOff),nOrb(iSym))
         End Do

*     Construct exchange contribution to Wij
*     --------------------------------------

         Do iVec1 = 1, NumVec
            iOff = nOrb(iSym)*nOccAll(iSym)*(iVec1-1)
            Call dGemm_('T','N',nOccAll(iSym),nOccAll(iSym),nOrb(iSym),
     &                 1.0d0,Wrk(kLip+iOff),nOrb(iSym),
     &                 Wrk(kVip+iOff),nOrb(iSym),1.0d0,
     &                 Wrk(kWij2),nOccAll(iSym))
         End Do
      End Do

      Call ChoMP2_OpenF(iClos,1,iSym)

*     Add SCF-density to MP2-density contribution
*     -------------------------------------------
      Do iSym1 = 1, nSym
         Do i = 1, nOcc(iSym1) + nFro(iSym1)
            Work(ipDensity(iSym1)+(i-1)*nOrb(iSym1)+(i-1)) =
     &        Work(ipDensity(iSym1)+(i-1)*nOrb(iSym1)+(i-1)) + 2.0d0
         End Do
      End Do

*     Construct Mp2 + SCF W-density
*     -----------------------------

      Do i = 1, nFro(iSym)
         iOrb = i
         Do j = 1, nFro(iSym)
            jOrb = j
            term = 0.0d0
            If(i .eq. j) term = Efro(iFro(iSym)+i)
            Work(ipWDensity(iSym) + jOrb-1 + (iOrb-1)*nOrb(iSym)) =
     &          - Wrk(kWJK(iSym) + j-1 + nFro(iSym)*(i-1)) + 2.0d0*term
            Work(ipWDensity(iSym) + jOrb-1 + (iOrb-1)*nOrb(iSym)) =
     &           Work(ipWDensity(iSym) + jOrb-1 + (iOrb-1)*nOrb(iSym)) -
     &           Wrk(kWij2 + j-1 + nOccAll(iSym)*(i-1))
         End Do
      End Do
      Do i = 1, nOcc(iSym)
            iOrb = nFro(iSym) + i
            Do j = 1, nFro(iSym)
               jOrb = j
               Work(ipWDensity(iSym) + jOrb-1 + (iOrb-1)*nOrb(iSym)) =
     &           - 0.5d0* Wrk(kWiK(iSym) + i-1 + nOcc(iSym)*(j-1))
               Work(ipWDensity(iSym) + jOrb-1 + (iOrb-1)*nOrb(iSym)) =
     &           Work(ipWDensity(iSym) + jOrb-1 + (iOrb-1)*nOrb(iSym))
     &         - Wrk(kWij2 + jOrb-1 + nOccAll(iSym)*(iOrb-1))
               Work(ipWDensity(iSym) + iOrb-1 + (jOrb-1)*nOrb(iSym)) =
     &           - 0.5d0* Wrk(kWiK(iSym) + i-1 + nOcc(iSym)*(j-1))
               Work(ipWDensity(iSym) + iOrb-1 + (jOrb-1)*nOrb(iSym)) =
     &           Work(ipWDensity(iSym) + iOrb-1 + (jOrb-1)*nOrb(iSym))
     &         - Wrk(kWij2 + jOrb-1 + nOccAll(iSym)*(iOrb-1))
            End Do
         End Do
      Do i = 1, nOcc(iSym)
         iOrb = nFro(iSym) + i
         Do j = 1, nOcc(iSym)
            jOrb = nFro(iSym) + j
            term = 0.0d0
            If(i .eq. j) term = EOcc(iOcc(iSym)+i)
            Work(ipWDensity(iSym) + jOrb-1 + (iOrb-1)*nOrb(iSym)) =
     &         - Wrk(kWij(iSym) + j-1 + nOcc(iSym)*(i-1)) + 2.0d0*term
            Work(ipWDensity(iSym) + jOrb-1 + (iOrb-1)*nOrb(iSym)) =
     &           Work(ipWDensity(iSym) + jOrb-1 + (iOrb-1)*nOrb(iSym))
     &         - Wrk(kWij2 + jOrb-1 + nOccAll(iSym)*(iOrb-1))
         End Do
      End Do

      Do i = 1, nFro(iSym)
         iOrb = i
         Do j = 1, nVir(iSym)
            jOrb = nFro(iSym) + nOcc(iSym) + j
            Work(ipWDensity(iSym) + jOrb-1 + (iOrb-1)*nOrb(iSym)) =
     &           -0.5d0*Wrk(kWaK(iSym) + j-1 + nVir(iSym)*(i-1))
            Work(ipWDensity(iSym) + iOrb-1 + (jOrb-1)*nOrb(iSym)) =
     &           -0.5d0*Wrk(kWaK(iSym) + j-1 + nVir(iSym)*(i-1))
         End Do
      End Do

      Do i = 1, nOcc(iSym)
         iOrb = nFro(iSym) + i
         Do j = 1, nVir(iSym)
            jOrb = nFro(iSym) + nOcc(iSym) + j
            Work(ipWDensity(iSym) + jOrb-1 + (iOrb-1)*nOrb(iSym)) =
     &           -0.5d0*Wrk(kWai(iSym) + j-1 + nVir(iSym)*(i-1))
            Work(ipWDensity(iSym) + iOrb-1 + (jOrb-1)*nOrb(iSym)) =
     &           -0.5d0*Wrk(kWai(iSym) + j-1 + nVir(iSym)*(i-1))
         End Do
      End Do


      Do i = 1, nVir(iSym)
         iOrb = nFro(iSym) + nOcc(iSym) +i
         Do j = 1, nVir(iSym)
            jOrb = nFro(iSym) + nOcc(iSym) + j
            Work(ipWDensity(iSym) + jOrb-1 + (iOrb-1)*nOrb(iSym)) =
     &         - Wrk(kWab(iSym) + j-1 + nVir(iSym)*(i-1))
         End Do
      End Do

#ifdef _DEBUGPRINT_
      Write(6,*) 'Full One-electron Density Matrix'
      iOff = 0
      Do iSym = 1, nSym
         Do i = 1, nOrb(iSym)*nOrb(iSym)
            Write(6,*) Work(ipMP2D + i-1+iOff)
         End Do
         iOff = iOff + nOrb(iSym)*nOrb(iSym)
      End Do
      Write(6,*) 'Full One-electron energy weighted D Matrix'
      iOff = 0
      Do iSym = 1, nSym
         Do i = 1, nOrb(iSym)*nOrb(iSym)
            Write(6,*) Work(ipMP2W + i-1+iOff)
         End Do
         iOff = iOff + nOrb(iSym)*nOrb(iSym)
      End Do
#endif


      Return
      End
