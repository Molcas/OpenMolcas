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
! Copyright (C) 2007, Thomas Bondo Pedersen                            *
!***********************************************************************
      SubRoutine ChoMP2_BackTra(iTyp,COcc,CVir,BaseName_AO,DoDiag,Diag)
!
!     Thomas Bondo Pedersen, Dec. 2007.
!
!     Purpose: Backtransform MO vectors to AO basis.
!              The result vectors are stored in lower triangular
!              storage.
!
!     Note: do not call this routine directly; use ChoMP2_VectorMO2AO()
!           instead !!
!
      use stdalloc
      Implicit None
      Integer iTyp
      Real*8  COcc(*), CVir(*)
      Character(LEN=3) BaseName_AO
      Logical DoDiag
      Real*8  Diag(*)

#include "cholesky.fh"
#include "choorb.fh"
#include "chomp2.fh"

      Character(LEN=14), Parameter:: SecNam = 'ChoMP2_BackTra'

      Character(LEN=4) FullName_AO

      Integer iSym, iSymb, iSyma, iSymi, iSymAl, iSymBe
      Integer nAB(8), iAB(8,8), nAB_Tot
      Integer lU_AO, l_Buf
      Integer MaxInCore, nVecInCore, nVecOnDisk, iVec, iOpt, iAdr, lVec
      Integer kCVir, kCOcc, kMOVec, kAOVec, kTemp, kDiag
      Integer na, ni, nAl, AlBe, kD

      Integer MulD2h, k, l

      Real*8, Allocatable:: AOVec(:), Temp(:), MOVec(:), Buf(:)

      MulD2h(k,l)=iEOr(k-1,l-1)+1

!     Set up index arrays.
!     --------------------

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

!     Backtransform.
!     --------------

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

         Call mma_allocate(AOVec,nAB(iSym),Label='AOVec')
         Call mma_allocate(Temp,nT1AOT(iSym),Label='Temp')
         Call mma_allocate(MOVec,nT1Am(iSym),Label='MOVec')

         Call mma_maxDBLE(l_Buf)
         If (l_Buf .lt. nAB(iSym)) Then
            Call ChoMP2_Quit(SecNam,'Insufficient memory!',' ')
         Else
            Call mma_allocate(Buf,l_Buf,Label='Buf')
         End If
         MaxInCore = min(l_Buf/nAB(iSym),nMP2Vec(iSym))

         nVecOnDisk = 0
         nVecInCore = 0
         Do iVec = 1,nMP2Vec(iSym)

            iOpt = 2
            iAdr = nT1Am(iSym)*(iVec-1) + 1
            lVec = nT1Am(iSym)
            Call ddaFile(lUnit_F(iSym,iTyp),iOpt,MOVec,lVec,            &
     &                   iAdr)

            Do iSymi = 1,nSym
               iSyma  = MulD2h(iSymi,iSym)
               iSymAl = iSyma
               kCVir  = iAOVir(iSymAl,iSyma) + 1
               kMOVec = 1 + iT1Am(iSyma,iSymi)
               kTemp  = 1 + iT1AOT(iSymi,iSymAl)
               na     = max(nVir(iSyma),1)
               nAl    = max(nBas(iSymAl),1)
               ni     = max(nOcc(iSymi),1)
               Call DGEMM_('T','T',nOcc(iSymi),nBas(iSymAl),nVir(iSyma),&
     &                    1.0d0,MOVec(kMOVec),na,CVir(kCVir),nAl,       &
     &                    0.0d0,Temp(kTemp),ni)
            End Do

            Do iSymBe = 1,nSym
               iSymAl = MulD2h(iSymBe,iSym)
               iSymi  = iSymBe
               kTemp  = 1 + iT1AOT(iSymi,iSymAl)
               kCOcc  = iT1AOT(iSymi,iSymBe) + 1
               kAOVec = 1 + iAB(iSymAl,iSymBe)
               ni     = max(nOcc(iSymi),1)
               nAl    = max(nBas(iSymAl),1)
               Call DGEMM_('T','N',                                     &
     &                    nBas(iSymAl),nBas(iSymBe),nOcc(iSymi),        &
     &                    1.0d0,Temp(kTemp),ni,COcc(kCOcc),ni,          &
     &                    0.0d0,AOVec(kAOVec),nAl)
            End Do

            If (DoDiag) Then
               Do AlBe = 1,nAB(iSym)
                  kD = kDiag + AlBe
                  Diag(kD) = Diag(kD) + AOVec(AlBe)**2
               End Do
            End If

            Call dCopy_(nAB(iSym),AOVec,1,Buf(1+nVecInCore),MaxInCore)
            nVecInCore = nVecInCore + 1

            If (nVecInCore.eq.MaxInCore .or.                            &
     &          iVec.eq.nMP2Vec(iSym)) Then
               Do AlBe = 1,nAB(iSym)
                  iOpt = 1
                  iAdr = nMP2Vec(iSym)*(AlBe-1) + nVecOnDisk + 1
                  lVec = nVecInCore
                  Call ddaFile(lU_AO,iOpt,                              &
     &                         Buf(1+MaxInCore*(AlBe-1)),lVec,iAdr)
               End Do
               nVecOnDisk = nVecOnDisk + nVecInCore
               nVecInCore = 0
            End If

         End Do
#if defined (_DEBUGPRINT_)
         If (nVecOnDisk.ne.nMP2Vec(iSym) .or. nVecInCore.ne.0) Then
            Call ChoMP2_Quit(SecNam,'Logical bug detected!',' [1]')
         End If
#endif

         Call mma_deallocate(Buf)
         Call mma_deallocate(MOVec)
         Call mma_deallocate(Temp)
         Call mma_deallocate(AOVec)

         Call daClos(lU_AO)
         iOpt = 2
         Call ChoMP2_OpenF(iOpt,iTyp,iSym)

  100    Continue ! cycle loop point

         If (DoDiag) kDiag = kDiag + nAB(iSym)

      End Do

      End
