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
      SubRoutine ChoMP2g_Reord_R(irc,Wrk,lWrk)
*
*      Jonas Bostrom, Apr 2010
*
*      Purpose: To reorder R-vectors so it is practical to access
*               one ia-piece at the time.

#include "implicit.fh"
#include "chomp2g.fh"
#include "chomp2.fh"
#include "chomp2_cfg.fh"
#include "cholesky.fh"
#include "choorb.fh"
#include "WrkSpc.fh"
*
      Character Fname*5
      Real*8 Wrk(lWrk)
*
      Character*7  ThisNm
      Character*15 SecNam
      Parameter (SecNam = 'ChoMP2g_Reord_r', ThisNm = 'Reord_r')
*
      MulD2h(i,j)=iEor(i-1,j-1) + 1
      iAdrVec(i,j,k) = (i-1) + (j-1)*nSym + (k-1)*nSym*nSym
*
      iTypR = 2
      iVecOV = 6
      maxvalue = 1000

*     Do not delete vectors
*     ---------------------
      iClos = 2

      iSeed = 7
      LuRInv(1) = IsFreeUnit(iSeed)
      Write(Fname,'(A4,I1)') 'TMPV',4
      Call DaName_MF_WA(LuRInv(1),Fname)
*
      LuRInv(2) = IsFreeUnit(iSeed)
      Write(Fname,'(A4,I1)') 'TMPV',5
      Call DaName_MF_WA(LuRInv(2),Fname)

      iAdr1 = 1
      iAdr2 = 1
      Do iSymI = 1, nSym
         Do iSymA = 1, nSym
            iSym = MulD2h(iSymA,iSymI)
            Do iI = 1, nOcc(iSymI)
               iWork(ipAdrR1 + iAdrVec(iSymA,iSymI,iI)) = iAdr1
               iAdr1 = iAdr1 + nVir(iSymA)*nMP2Vec(iSym)
            End Do
            Do iA = 1, nVir(iSymA)
               iWork(ipAdrR2 + iAdrVec(iSymA,iSymI,iA)) = iAdr2
               iAdr2 = iAdr2 + nOcc(iSymI)*nMP2Vec(iSym)
            End Do
         End Do
      End Do
*
      Do iSym = 1, nSym
         If(nMP2Vec(iSym) .eq. 0) Go To 10
         nVec = Min(maxvalue,nMP2Vec(iSym))
         If(nVec .lt. 1) Then
            Call ChoMP2_Quit(SecNam,'Insufficient memory','[1]')
         End If
         nBatR = (nMP2Vec(iSym)-1)/nVec + 1

*        Allocate memory for Ria-vectors
*        -------------------------------

         lRia = nMoMo(iSym,iVecOV)*nVec
         kRia1 = 1
         kEndRia1 = kRia1 + lRia

         kRia2 = kEndRia1
         kEndRia2 = kRia2 + lRia

         kRia3 = kEndRia2

*        Open Cholesky amplitude vectors
*        -------------------------------
         Call ChoMP2_OpenF(1,iTypR,iSym)

         Do iBat = 1, nBatR
            If(iBat .eq. nBatR) Then
               NumVec = nMP2Vec(iSym) - nVec*(nBatR-1)
            Else
               NumVec = nVec
            End If
            iVec = nVec*(iBat-1) + 1

*           Read Amplitude vectors
*           ----------------------
            iOpt = 2
            lTot = nMoMo(iSym,iVecOV)*NumVec
            iAdr = nMoMo(iSym,iVecOV)*(iVec-1) + 1
            Call dDaFile(lUnit_F(iSym,iTypR),iOpt,Wrk(kRia1),lTot,iAdr)

            Do iVec1 = 1, NumVec
               Do iSymI = 1, nSym
                  iSymA = MulD2h(iSymI,iSym)
                  Do iI = 1, nOcc(iSymI)
                     ioffset1 = (iI-1)*nVir(iSymA)+iT1am(iSymA,iSymI) +
     &                         (iVec1-1)*nMoMo(iSym,iVecOV)
                     iOffset2 = (iVec1-1)*nVir(iSymA) +
     &                          (iI-1)*NumVec*nVir(iSymA) +
     &                           iT1am(iSymA,iSymI)*NumVec
                     Call dCopy_(nVir(iSymA),Wrk(kRia1+iOffset1),1,
     &                                      Wrk(kRia2+iOffset2),1)
                  End Do
               End Do
            End Do

*           Put the reordered vectors on disk
            iOpt = 1
            Do iSymI = 1, nSym
               iSymA = MulD2h(iSymI,iSym)
               Do iI = 1, nOcc(iSymI)
                  lTot = nVir(iSymA)*NumVec
                  iAdr = iWork(ipAdrR1 + iAdrVec(iSymA,iSymI,iI))
     &                 + (iVec-1)*nVir(iSymA)
                  iOffset2 = (iI-1)*NumVec*nVir(iSymA) +
     &                        iT1am(iSymA,iSymI)*NumVec
                  Call dDaFile(LuRInv(1),iOpt, Wrk(kRia2+iOffSet2),
     &                         lTot,iAdr)
               End Do
            End Do



         End Do                 !iBat
         Do iBat = 1, nBatR
            If(iBat .eq. nBatR) Then
               NumVec = nMP2Vec(iSym) - nVec*(nBatR-1)
            Else
               NumVec = nVec
            End If
            iVec = nVec*(iBat-1) + 1

*           Read Amplitude vectors
*           ----------------------
            iOpt = 2
            lTot = nMoMo(iSym,iVecOV)*NumVec
            iAdr = nMoMo(iSym,iVecOV)*(iVec-1) + 1
            Call dDaFile(lUnit_F(iSym,iTypR),iOpt,Wrk(kRia1),lTot,iAdr)
*
            Do iSymI = 1, nSym
               iSymA = MulD2h(iSymI,iSym)
               Do iI = 1, nOcc(iSymI)
                  Do iA = 1, nVir(iSymA)

                     iOffSet1 = iA-1 + (iI-1)*nVir(iSymA)
     &                        + iT1am(iSymA,iSymI)
                     iOffSet2 = iI-1 + (iA-1)*NumVec*nOcc(iSymI)
     &                        + iT1am(iSymA,iSymI)*NumVec
                     Call dCopy_(NumVec,Wrk(kRia1+iOffset1),
     &                          nMoMo(iSym,iVecOV),Wrk(kRia2+iOffset2),
     &                          nOcc(iSymI))
                  End Do
               End Do
            End Do
*
*           Put the reordered vectors on disk
            iOpt = 1
            Do iSymI = 1, nSym
               iSymA = MulD2h(iSymI,iSym)
               Do iA = 1, nVir(iSymA)
                  lTot = nOcc(iSymI)*NumVec
                  iAdr = iWork(ipAdrR2 + iAdrVec(iSymA,iSymI,iA))
     &                 + (iVec-1)*nOcc(iSymI)
                  iOffset2 = (iA-1)*NumVec*nOcc(iSymI) +
     &                        iT1am(iSymA,iSymI)*NumVec
                  Call dDaFile(LuRInv(2),iOpt, Wrk(kRia2+iOffset2),
     &                         lTot,iAdr)
               End Do
            End Do



         End Do                 !iBat

         Call ChoMP2_OpenF(iClos,iTypR,iSym)


 10      Continue
      End Do !iSym

      Call DaClos(LuRInv(1))
      Call DaClos(LuRInv(2))

c Avoid unused argument warnings
      If (.False.) Call Unused_integer(irc)
      End
