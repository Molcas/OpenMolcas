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
* Copyright (C) 2007, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine ChoMP2_BackTra(iTyp,COcc,CVir,BaseName_AO,DoDiag,Diag)
C
C     Thomas Bondo Pedersen, Dec. 2007.
C
C     Purpose: Backtransform MO vectors to AO basis.
C              The result vectors are stored in lower triangular
C              storage.
C
C     Note: do not call this routine directly; use ChoMP2_VectorMO2AO()
C           instead !!
C
      Implicit None
      Integer iTyp
      Real*8  COcc(*), CVir(*)
      Character*3 BaseName_AO
      Logical DoDiag
      Real*8  Diag(*)
#include "cholesky.fh"
#include "choorb.fh"
#include "chomp2.fh"
#include "WrkSpc.fh"

      Character*14 SecNam
      Parameter (SecNam = 'ChoMP2_BackTra')

      Character*4 FullName_AO

      Integer iSym, iSymb, iSyma, iSymi, iSymAl, iSymBe
      Integer nAB(8), iAB(8,8), nAB_Tot
      Integer lU_AO
      Integer ip_AOVec, l_AOVec, ip_Temp, l_Temp, ip_MOVec, l_MOVec
      Integer ip_Buf, l_Buf
      Integer MaxInCore, nVecInCore, nVecOnDisk, iVec, iOpt, iAdr, lVec
      Integer kCVir, kCOcc, kMOVec, kAOVec, kAOV, kTemp, kDiag
      Integer na, ni, nAl, AlBe, kD, kA

      Integer MulD2h, k, l
      MulD2h(k,l)=iEOr(k-1,l-1)+1

C     Set up index arrays.
C     --------------------

      Call iCopy(8*8,[0],0,iAB,1)
      nAB_Tot = 0
      Do iSym = 1,nSym
         nAB(iSym) = 0
         Do iSymb = 1,nSym
            iSyma = MulD2h(iSymb,iSym)
            iAB(iSyma,iSymb) = nAB(iSym)
            nAB(iSym) = nAB(iSym) + nBas(iSyma)*nBas(iSymb)
         End Do
         nAB_Tot = nAB_Tot + nAB(iSym)
      End Do

C     Backtransform.
C     --------------

      If (DoDiag) Then
         Call dCopy_(nAB_Tot,[0.0d0],0,Diag,1)
      End If

      kDiag = 0
      Do iSym = 1,nSym

         If (nAB(iSym) .lt. 1) Go To 100     ! cycle loop
         If (nMP2Vec(iSym) .lt. 1) Go To 100 ! cycle loop

         iOpt = 1
         Call ChoMP2_OpenF(iOpt,iTyp,iSym)
         Write(FullName_AO,'(A3,I1)') BaseName_AO,iSym
         lU_AO = 7
         Call daName_MF_WA(lU_AO,FullName_AO)

         l_AOVec = nAB(iSym)
         l_Temp  = nT1AOT(iSym)
         l_MOVec = nT1Am(iSym)
         Call GetMem('AOVec','Allo','Real',ip_AOVec,l_AOVec)
         Call GetMem('Temp','Allo','Real',ip_Temp,l_Temp)
         Call GetMem('MOVec','Allo','Real',ip_MOVec,l_MOVec)

         Call GetMem('GetMx','Max ','Real',ip_Buf,l_Buf)
         If (l_Buf .lt. nAB(iSym)) Then
            Call ChoMP2_Quit(SecNam,'Insufficient memory!',' ')
         Else
            Call GetMem('Buffer','Allo','Real',ip_Buf,l_Buf)
         End If
         MaxInCore = min(l_Buf/nAB(iSym),nMP2Vec(iSym))

         nVecOnDisk = 0
         nVecInCore = 0
         Do iVec = 1,nMP2Vec(iSym)

            iOpt = 2
            iAdr = nT1Am(iSym)*(iVec-1) + 1
            lVec = nT1Am(iSym)
            Call ddaFile(lUnit_F(iSym,iTyp),iOpt,Work(ip_MOVec),lVec,
     &                   iAdr)

            Do iSymi = 1,nSym
               iSyma  = MulD2h(iSymi,iSym)
               iSymAl = iSyma
               kCVir  = iAOVir(iSymAl,iSyma) + 1
               kMOVec = ip_MOVec + iT1Am(iSyma,iSymi)
               kTemp  = ip_Temp + iT1AOT(iSymi,iSymAl)
               na     = max(nVir(iSyma),1)
               nAl    = max(nBas(iSymAl),1)
               ni     = max(nOcc(iSymi),1)
               Call DGEMM_('T','T',nOcc(iSymi),nBas(iSymAl),nVir(iSyma),
     &                    1.0d0,Work(kMOVec),na,CVir(kCVir),nAl,
     &                    0.0d0,Work(kTemp),ni)
            End Do

            Do iSymBe = 1,nSym
               iSymAl = MulD2h(iSymBe,iSym)
               iSymi  = iSymBe
               kTemp  = ip_Temp + iT1AOT(iSymi,iSymAl)
               kCOcc  = iT1AOT(iSymi,iSymBe) + 1
               kAOVec = ip_AOVec + iAB(iSymAl,iSymBe)
               ni     = max(nOcc(iSymi),1)
               nAl    = max(nBas(iSymAl),1)
               Call DGEMM_('T','N',
     &                    nBas(iSymAl),nBas(iSymBe),nOcc(iSymi),
     &                    1.0d0,Work(kTemp),ni,COcc(kCOcc),ni,
     &                    0.0d0,Work(kAOVec),nAl)
            End Do

            If (DoDiag) Then
               kAOV = ip_AOVec - 1
               Do AlBe = 1,nAB(iSym)
                  kD = kDiag + AlBe
                  kA = kAOV + AlBe
                  Diag(kD) = Diag(kD) + Work(kA)*Work(kA)
               End Do
            End If

            Call dCopy_(nAB(iSym),Work(ip_AOVec),1,
     &                           Work(ip_Buf+nVecInCore),MaxInCore)
            nVecInCore = nVecInCore + 1

            If (nVecInCore.eq.MaxInCore .or.
     &          iVec.eq.nMP2Vec(iSym)) Then
               Do AlBe = 1,nAB(iSym)
                  iOpt = 1
                  iAdr = nMP2Vec(iSym)*(AlBe-1) + nVecOnDisk + 1
                  lVec = nVecInCore
                  Call ddaFile(lU_AO,iOpt,
     &                         Work(ip_Buf+MaxInCore*(AlBe-1)),lVec,
     &                         iAdr)
               End Do
               nVecOnDisk = nVecOnDisk + nVecInCore
               nVecInCore = 0
            End If

         End Do
#if defined (_DEBUG_)
         If (nVecOnDisk.ne.nMP2Vec(iSym) .or. nVecInCore.ne.0) Then
            Call ChoMP2_Quit(SecNam,'Logical bug detected!',' [1]')
         End If
#endif

         Call GetMem('Buffer','Free','Real',ip_Buf,l_Buf)
         Call GetMem('MOVec','Free','Real',ip_MOVec,l_MOVec)
         Call GetMem('Temp','Free','Real',ip_Temp,l_Temp)
         Call GetMem('AOVec','Free','Real',ip_AOVec,l_AOVec)

         Call daClos(lU_AO)
         iOpt = 2
         Call ChoMP2_OpenF(iOpt,iTyp,iSym)

  100    Continue ! cycle loop point

         If (DoDiag) kDiag = kDiag + nAB(iSym)

      End Do

      End
