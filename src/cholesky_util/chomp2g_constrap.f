************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************

      Subroutine ChoMP2g_ConstrAP(irc,Scr,lScr,Type,
     &                            iSym,nVec,Ap,lAp,Dens,lDens,factor)


#include "implicit.fh"
#include "chomp2g.fh"
#include "chomp2.fh"
#include "chomp2_cfg.fh"
#include "cholesky.fh"
#include "choorb.fh"
#include "WrkSpc.fh"
*
      Real*8 Scr(lScr), Ap(lAp), Dens(lDens)

      Character  Type(1:4)
      Logical DoX, DoY, DoZ
      Integer iOffp(8), iOffAp(8)
      Integer nOrb1(8,4), iOrbType(4)
*
      Character*8  ThisNm
      Character*16 SecNam
      Parameter (SecNam = 'ChoMP2g_ConstrAP', ThisNm = 'ConstrAP')
*************************************
      MulD2h(i,j)=iEor(i-1,j-1) + 1
*************************************

*     Setup some lengths
*     ------------------
      nTypes = 4

      iTypL = 1

*     Setup lengths of needed cholesky/intermediate vectors
*     ----------------------------------------------------

*     Read the char-array Type = 'pqri'
*     ---------------------------------
      Do i = 1, nTypes
         If(Type(i) .eq. 'f') Then
            iOrbType(i) = 1
            Do iSym1 = 1, nSym
               nOrb1(iSym1,i) = nFro(iSym1)
            End Do
         Else If (Type(i) .eq. 'o') Then
            iOrbType(i) = 2
            Do iSym1 = 1, nSym
               nOrb1(iSym1,i) = nOcc(iSym1)
            End Do
         Else If (Type(i) .eq. 'v') Then
            iOrbType(i) = 3
            Do iSym1 = 1, nSym
               nOrb1(iSym1,i) = nVir(iSym1)
            End Do
         Else
            Write(6,*) 'Forbidden Type pqrs in', SecNam
            irc = -1
            Return
         End If
      End Do

*     Set vector type index
*     ---------------------
      iPQ = 3*(iOrbType(1)-1) + iOrbType(2)
      iIR = 3*(iOrbType(4)-1) + iOrbType(3)
      iIQ = 3*(iOrbType(4)-1) + iOrbType(2)
      iRP = 3*(iOrbType(3)-1) + iOrbType(1)
      iIP = 3*(iOrbType(4)-1) + iOrbType(1)
      iRQ = 3*(iOrbType(3)-1) + iOrbType(2)

*     Set vector lengths
*     ------------------
      nPQ = nMoMo(iSym,iPQ)
      nIR = nMoMo(iSym,iIR)
      nIQ = nMoMo(iSym,iIQ)
      nRP = nMoMo(iSym,iRP)
      nIP = nMoMo(iSym,iIP)
      nRQ = nMoMo(iSym,iRQ)

*     Set adress-offset for cholesky vector adress
*     --------------------------------------------
      iAdrLpq = iAdrOff(iSym,iPQ)
      iAdrLir = iAdrOff(iSym,iIR)
      iAdrLiq = iAdrOff(iSym,iIQ)
      iAdrLrp = iAdrOff(iSym,iRP)
      iAdrLip = iAdrOff(iSym,iIP)
      iAdrLrq = iAdrOff(iSym,iRQ)

*     Set
*

#ifdef _DEBUGPRINT_
      Write(6,*) 'iPQ', iPQ
      Write(6,*) 'iIQ', iIQ
      Write(6,*) 'iIR', iIR
      Write(6,*) 'iIP', iIP
      Write(6,*) 'iRQ', iRQ
      Write(6,*) 'iRP', iRP
*
      Write(6,*) 'nPQ', nPQ
      Write(6,*) 'nIQ', nIQ
      Write(6,*) 'nIR', nIR
      Write(6,*) 'nIP', nIP
      Write(6,*) 'nRQ', nRQ
      Write(6,*) 'nRP', nRP
#endif
*     Decide what type of intermediate vectors is needed
*     --------------------------------------------------
      DoX = .True.
      DoY = .True.
      DoZ = .True.
*
      If(iSym.ne.1) DoX = .False.
      If((nPQ .eq. 0) .or. (nIR .eq. 0)) DoX = .False.
      If((nRP .eq. 0) .or. (nIQ .eq. 0)) DoY = .False.
      If((nIP .eq. 0) .or. (nRQ .eq. 0)) DoZ = .False.
      If(Type(1) .eq. Type(2)) DoZ = .False.

      If(DoY.or.DoZ) Then
         iOffp(1) = 0
         iOffap(1) = 0
         Do iSym1 = 2, nSym
            iOffp(iSym1) = iOffp(iSym1-1) +
     &                      nOrb1(iSym1-1,1)*nOrb1(iSym1-1,2)
            iOffap(iSym1) = iOffap(iSym1-1) +
     &                      nOrb1(iSym1-1,3)*nOrb1(iSym1-1,4)
         End Do
      End If

*     Allocate memory for L-vectors
*     ------------------------------

*        Lri:
         lLir = nIR*nVec
         kLir = 1
         kEndLir = kLir +lLir

*        Lpq:
         lLpq = nPQ*nVec
         kLpq = kEndLir
         kEndLpq = kLpq + lLpq

*        Lrp:
         lLrp = nRP*nVec
         kLrp = kEndLpq
         kEndLrp = kLrp + lLrp

*        Lrq:
         lLrq = nRQ*nVec
         kLrq = kEndLrp
         kEndLrq = kLrq + lLrq

*        Lip:
         lLip = nIP*nVec
         kLip = kEndLrq
         kEndLip = kLip + lLip

*        Liq:
         lLiq = nIQ*nVec
         kLiq = kEndLip
         kEndLiq = kLiq + lLiq

*     Allocate memory for Intermediate Vectors
*     ----------------------------------------
*        X-vector
         lX = NumCho(iSym)
         kX = kEndLiq
         kEndX = kX + lX
*        Y-vector
         lY = nVec*nIP
         kY = kEndX
         kEndY = kY + lY
*        Z-vector
         lZ = nVec*nIQ
         kZ = kEndY
         kEndZ = kZ + lZ
*
      nBatL = (NumCho(iSym)-1)/nVec + 1

*     Construction of Intermediate vectors
*     ------------------------------------
      Do iBat = 1, nBatL
         If (iBat .eq. nBatL) Then
            NumVec = NumCho(iSym) - nVec*(nBatL-1)
         Else
            NumVec = nVec
         End If
         iVec = nVec*(iBat-1) + 1

*        Read Lpq-vectors
*        ----------------
         If(DoX) Then
            iOpt = 2
            lTot = nPQ*NumVec
            iAdr = nPQ*(iVec-1) + 1 +
     &             iAdrLpq
            Call dDaFile(lUnit_F(iSym,iTypL),iOpt,Scr(kLpq),
     &                   lTot,iAdr)
         End If

*        Read Liq-vectors
*        ----------------
         If(DoY) Then
            iOpt = 2
            lTot = nIQ*NumVec
            iAdr = nIQ*(iVec-1) + 1 +
     &           iAdrLiq
            Call dDaFile(lUnit_F(iSym,iTypL),iOpt,Scr(kLiq),
     &                   lTot,iAdr)
         End If

*
         If(DoZ) Then
            iOpt = 2
            lTot = nIP*NumVec
            iAdr = nIP*(iVec-1) + 1 + iAdrLip
            Call dDaFile(lUnit_F(iSym,iTypL),iOpt,Scr(kLip),
     &                   lTot,iAdr)
         End If

*        Construct X-vector
*        ------------------
         If(DoX) Then
            nRow = 1
            Call dGemm_('T','N',nRow,NumVec,nPQ,
     &                 1.0d0,Dens(1),nPQ,
     &                       Scr(kLpq),nPQ,0.0d0,
     &                       Scr(kX+iVec-1),nRow)
         End If

*        Construct Y-vector
*        ------------------
         If(DoY) Then
            Do iVec1 = 1, NumVec
               iOffL1 = 0
               iOffY1 = 0
               Do iSymI = 1, nSym
                  iSymP = MulD2h(iSym,iSymI)
                  iSymQ = iSymP
                  nI = nOrb1(iSymI,4)
                  nP = nOrb1(iSymP,1)
                  nQ = nOrb1(iSymQ,2)
                  nR = nOrb1(iSymI,3)
                  If(nQ*nP*nI*nR .eq. 0) Go To 101
                  iOffL = (iVec1-1)*nIQ + iOffL1
                  iOffY = (iVec1-1)*nIP + iOffY1
                  Call dGemm_('T','N', nP, nI, nQ,1.0d0,
     &                       Dens(1+iOffP(iSymP)), nQ,
     &                       Scr(kLiq+iOffL), nQ, 0.0d0,
     &                       Scr(kY+iOffY), nP)
 101              Continue
                  iOffL1 = iOffL1 + nQ*nI
                  iOffY1 = iOffY1 + nP*nI
               End Do
            End Do
            If(nBatL .ne. 1) Then
               iOpt = 1
               lTot = nIP*NumVec
               iAdr = nIP*(iVec-1)+1
               Call dDaFile(LuVVec, iOpt, Scr(kY),lTot,iAdr)
            End If
         End If
*        Construct Z-vectors
*        -------------------
         If(DoZ) Then
            Do iVec1 = 1, NumVec
               iOffL1 = 0
               iOffZ1 = 0
               Do iSymI = 1, nSym
                  iSymP = MulD2h(iSym,iSymI)
                  iSymQ = iSymP
                  nI = nOrb1(iSymI,4)
                  nP = nOrb1(iSymP,1)
                  nQ = nOrb1(iSymQ,2)
                  nR = nOrb1(iSymI,3)
                  If(nQ*nP*nI*nR .eq. 0) Go To 102
                  iOffL = (iVec1-1)*nIP + iOffL1
                  iOffZ = (iVec1-1)*nIQ + iOffZ1
                  Call dGemm_('N','N', nQ, nI, nP,1.0d0,
     &                       Dens(1+iOffP(iSymQ)), nQ,
     &                       Scr(kLip+iOffL), nP, 0.0d0,
     &                       Scr(kZ+iOffZ), nQ)
 102              Continue
                  iOffL1 = iOffL1 + nP*nI
                  iOffZ1 = iOffZ1 + nQ*nI
               End Do
            End Do
            If(nBatL .ne. 1) Then
               iOpt = 1
               lTot = nIQ*NumVec
               iAdr = nIQ*(iVec-1)+1
               Call dDaFile(LuWVec, iOpt, Scr(kZ),lTot,iAdr)
            End If
*
         End If
      End Do

*     Contracting to produce the final contribution to A*p
*     ----------------------------------------------------
      Do iBat = 1, nBatL
         If (iBat .eq. nBatL) Then
            NumVec = NumCho(iSym) - nVec*(nBatL-1)
         Else
            NumVec = nVec
         End If
         iVec = nVec*(iBat-1) + 1

*        Read Lir-vectors
*        ----------------
         If(DoX) Then
            iOpt = 2
            lTot = nIR*NumVec
            iAdr = nIR*(iVec-1) + 1 +
     &             iAdrLir
            Call dDaFile(lUnit_F(iSym,iTypL),iOpt,Scr(kLir),
     &                   lTot,iAdr)
         End If

*        Read Lrp-vectors
*        ----------------
         If(DoY) Then
            iOpt = 2
            lTot = nRP*NumVec
            iAdr = nRP*(iVec-1) + 1 + iAdrLrp
            Call dDaFile(lUnit_F(iSym,iTypL),iOpt,Scr(kLrp),
     &                   lTot,iAdr)
         End If

*        Read Lrq-vectors
*        ----------------
         If(DoZ) Then
            iOpt = 2
            lTot = nRQ*NumVec
            iAdr = nRQ*(iVec-1) + 1 + iAdrLrq
            Call dDaFile(lUnit_F(iSym,iTypL),iOpt,Scr(kLrq),
     &                   lTot,iAdr)
         End If
*
         If(DoX) Then
            iCol = 1
            Call dGemm_('N','N',nIR,iCol,NumVec,
     &                 4.0d0*factor, Scr(kLir),nIR,
     &                 Scr(kX+iVec-1),NumVec, 1.0d0,
     &                 Ap(1),nIR)
         End If
*
         If(DoY) Then
*           Read Y-vector (if they are stored on disk and not in memory)
*           ------------------------------------------------------------
            If(nBatL .ne. 1) Then
               iOpt = 2
               lTot = nIP*NumVec
               iAdr = nIP*(iVec-1)+1
               Call dDaFile(LuVVec, iOpt, Scr(kY),lTot,iAdr)
            End If
*
            yfactor = 1.0d0*factor
            If(.not. DoZ) yfactor = yfactor*2.0d0
            Do iVec1 = 1, NumVec
                  iOffL1 = 0
                  iOffY1 = 0
                  Do iSymI = 1, nSym
                     iSymP = MulD2h(iSym,iSymI)
                     iSymR = iSymI
                     nI = nOrb1(iSymI,4)
                     nR = nOrb1(iSymR,3)
                     nP = nOrb1(iSymP,1)
                     nQ = nOrb1(iSymP,2)
*
                     If(nI*nR*nP*nQ .eq. 0) Go To 201
                     iOffL = (iVec1-1)*nRP + iOffL1
                     iOffY = (iVec1-1)*nIP + iOffY1
                     Call dGemm_('T','N',nR, nI, nP,-1.0d0*yfactor,
     &                          Scr(kLrp+iOffL), nP,
     &                          Scr(kY+iOffY),nP ,1.0d0,
     &                          Ap(1+ iOffAP(iSymI)),nR)
 201                 Continue
                     iOffL1 = iOffL1 + nR*nP
                     iOffY1 = iOffY1 + nI*nP
                  End Do
               End Do


         End If

         If(DoZ) Then
*           Read Z-vector (if they are stored on disk and not in memory)
*           ------------------------------------------------------------
            If(nBatL .ne. 1) Then
               iOpt = 2
               lTot = nIQ*NumVec
               iAdr = nIQ*(iVec-1)+1
               Call dDaFile(LuWVec, iOpt, Scr(kZ),lTot,iAdr)
            End If
*
            Do iVec1 = 1, NumVec
                  iOffL1 = 0
                  iOffZ1 = 0
                  Do iSymI = 1, nSym
                     iSymQ = MulD2h(iSym,iSymI)
                     iSymR = iSymI
                     nI = nOrb1(iSymI,4)
                     nR = nOrb1(iSymR,3)
                     nQ = nOrb1(iSymQ,2)
                     nP = nOrb1(iSymQ,1)
                     If(nI*nR*nQ*nP .eq. 0) Go To 202
                     iOffL = (iVec1-1)*nRQ + iOffL1
                     iOffZ = (iVec1-1)*nIQ + iOffZ1
                     Call dGemm_('T','N',nR, nI, nQ,-1.0d0*factor,
     &                          Scr(kLrq+iOffL), nQ,
     &                          Scr(kZ+iOffZ),nQ ,1.0d0,
     &                          Ap(1+ iOffAP(iSymI)),nR)
 202                 Continue
                     iOffL1 = iOffL1 + nR*nQ
                     iOffZ1 = iOffZ1 + nI*nQ
                  End Do
               End Do
         End If
      End Do
      Return
      End
