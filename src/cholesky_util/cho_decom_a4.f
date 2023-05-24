!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      SubRoutine Cho_Decom_A4(Diag,LstQSP,NumSP,iPass)
!
!     Purpose: decompose qualified columns ("parallel" algorithm).
!
      use ChoArr, only: LQ_Tot, LQ
      use ChoVecBuf, only: nVec_in_Buf
      use stdalloc

      Implicit Real*8 (a-h,o-z)
      Real*8  Diag(*)
      Integer LstQSP(NumSP)
#include "cholesky.fh"
#include "choprint.fh"

      Character(LEN=12), Parameter:: SecNam = 'Cho_Decom_A4'

      Integer NumCho_Old(8), nQual_Old(8)
      Integer NumV(8), nkVec(8)

      Real*8, Allocatable:: KVScr(:), MQ(:), KVec(:), QDiag(:), xInt(:),
     &                      Wrk1(:)
      Integer, Allocatable:: IDKVec(:), iQScr(:)
!                                                                      *
!***********************************************************************
!                                                                      *
      Interface
      SubRoutine Cho_P_GetLQ(QVec,l_QVec,LstQSP,nQSP)
      Integer l_QVec, nQSP
      Real*8, Target::  QVec(l_Qvec)
      Integer LstQSP(nQSP)
      End SubRoutine Cho_P_GetLQ
      End Interface
!                                                                      *
!***********************************************************************
!                                                                      *

!     Print header.
!     -------------

      LenLin = 0
      If (iPrint .ge. Inf_Progress) Then
         Call Cho_Head(SecNam//
     &                 ': Decomposition of Qualified Diagonals','=',
     &                 80,LUPRI)
         Write(Lupri,'(/,A,I5,A,I4,A)')
     &   'Integral pass number',iPass,' (',NumSP,
     &   ' shell pair distributions calculated)'
         Write(Lupri,'(A,8I8)')
     &   '#Cholesky vec.: ',(NumCho(iSym),iSym=1,nSym)
         Write(Lupri,'(A,8I8)')
     &   '#vec. in buff.: ',(nVec_in_Buf(iSym),iSym=1,nSym)
         Write(Lupri,'(A,8I8)')
     &   '#qualified    : ',(nQual(iSym),iSym=1,nSym)
         Write(Lupri,'(A,8I8)')
     &   'Current  dim. : ',(nnBstr(iSym,2),iSym=1,nSym)
         Write(Lupri,'(A,8I8)')
     &   'Original dim. : ',(nnBstr(iSym,1),iSym=1,nSym)
         Write(Lupri,'(/,A,/,A,A)')
     &   '           #Vectors             Treated Diagonal',
     &   'Sym.     Sym.     Total     Index     Before      After',
     &   '   Conv. Neg.   New Max'
         LenLin = 79
         Write(Lupri,'(80A)') ('-',I=1,LenLin)
         Call Cho_Flush(Lupri)
         Call iCopy(nSym,NumCho,1,NumCho_Old,1)
      Else If (iPrint .ge. Inf_Pass) Then
         Write(Lupri,'(/,A,I4)')
     &   'Number of shell pair distributions calculated:',NumSP
         Write(Lupri,'(A,8I8)')
     &   '#Cholesky vec.: ',(NumCho(iSym),iSym=1,nSym)
         Write(Lupri,'(A,8I8)')
     &   '#vec. in buff.: ',(nVec_in_Buf(iSym),iSym=1,nSym)
         Write(Lupri,'(A,8I8)')
     &   '#qualified    : ',(nQual(iSym),iSym=1,nSym)
         Call Cho_Flush(Lupri)
         Call iCopy(nSym,NumCho,1,NumCho_Old,1)
      End If

!     Allocations.
!     ------------

      Call Cho_P_GetGV(numV,nSym)

      l_KVec = nQual(1)**2
      l_IDKVec = nQual(1)
      l_LQ = nQual(1)*NumV(1)
      Do iSym = 2,nSym
         l_KVec = l_KVec + nQual(iSym)**2
         l_IDKVec = l_IDKVec + nQual(iSym)
         l_LQ = l_LQ + nQual(iSym)*NumV(iSym)
      End Do
      l_LQ = max(l_LQ,1) ! because there might not be any prev. vecs.
      Call mma_allocate(KVec,l_KVec,Label='KVec')
      Call mma_allocate(IDKVec,l_IDKVec,Label='IDKVec')
      Call mma_allocate(QDiag,l_IDKVec,Label='QDiag')

!     Extract elements corresponding to qualified diagonals from
!     previous Cholesky vectors (if any).
!     ----------------------------------------------------------

      Call Cho_Timer(C1,W1)
      Call mma_allocate(LQ_Tot,l_LQ,Label='LQ_Tot')

      iEn = 0
      iSt = 1
      Do iSym = 1,nSym
         If (nQual(iSym)*NumV(iSym)>0) Then
            iEn = iEn + nQual(iSym)*NumV(iSym)
            LQ(iSym)%Array(1:nQual(iSym),1:NumV(iSym)) =>
     &                LQ_Tot(iSt:iEn)
            iSt = iEn + 1
         Else
            LQ(iSym)%Array => Null()
         End If
      End Do

      Call Cho_P_GetLQ(LQ_Tot,l_LQ,LstQSP,NumSP)
      Call Cho_Timer(C2,W2)
      tDecom(1,2) = tDecom(1,2) + C2 - C1
      tDecom(2,2) = tDecom(2,2) + W2 - W1

!     Extract qualified diagonal integral block.
!     ------------------------------------------

      Call Cho_Timer(C1,W1)

      Call mma_allocate(MQ,l_KVec,Label='MQ')
      Call Cho_P_GetMQ(MQ,SIZE(MQ),LstQSP,NumSP)

!     Decompose qualified diagonal block.
!     The qualified diagonals are returned in QDiag.
!     ----------------------------------------------

      Call Cho_Dec_Qual(Diag,LQ_Tot,MQ,KVec,IDKVec,nKVec,QDiag)

!     Deallocate MQ.
!     --------------

      Call mma_deallocate(MQ)

!     Reorder the elements of the K-vectors according to IDK ordering.
!     ----------------------------------------------------------------

      MxQ = nQual(1)
      Do iSym = 2,nSym
         MxQ = max(MxQ,nQual(iSym))
      End Do

      Call mma_allocate(KVScr,MxQ,Label='KVScr')

      kK1 = 0
      kK2 = kK1
      kID = 0
      Do iSym = 1,nSym
         Do iK = 1,nKVec(iSym)
            kK_1 = kK1 + nQual(iSym)*(iK-1) + 1
            Call dCopy_(nQual(iSym),KVec(kK_1),1,KVScr,1)
            kK_2 = kK2 + nKVec(iSym)*(iK-1)
            Do jK = 1,nKVec(iSym)
               lK = IDKVec(kID+jK)
               KVec(kK_2+jK) = KVScr(lK)
            End Do
         End Do
         kK1 = kK1 + nQual(iSym)**2
         kK2 = kK2 + nKVec(iSym)**2
         kID = kID + nQual(iSym)
      End Do

!     Reorder QDiag to IDK ordering.
!     ------------------------------

      kID = 0
      kQD = 0
      Do iSym = 1,nSym
         Call dCopy_(nQual(iSym),QDiag(kQD+1),1,KVScr,1)
         Do iK = 1,nKVec(iSym)
            lK = IDKVec(kID+iK)
            QDiag(kQD+iK) = KVScr(lK)
         End Do
         kQD = kQD + nQual(iSym)
         kID = kID + nQual(iSym)
      End Do

!     Reorder elements of LQ vectors to IDK ordering.
!     -----------------------------------------------

      kID = 0
      Do iSym = 1,nSym
         If (nQual(iSym)<1) Cycle
         Do jVec = 1,NumV(iSym)
            Call dCopy_(nQual(iSym),LQ(iSym)%Array(:,jVec),1,
     &                              KVScr,1)
            Do iK = 1,nKVec(iSym)
               lK = IDKVec(kID+iK)
               LQ(iSym)%Array(iK,jVec) = KVScr(lK)
            End Do
         End Do
         kID = kID + nQual(iSym)
      End Do

      Call mma_deallocate(KVScr)

!     Reset qualification index arrays to IDK ordering.
!     Local as well as global are reordered.
!     -------------------------------------------------

      Call iCopy(nSym,nQual,1,nQual_Old,1)
      Call mma_allocate(iQScr,MxQ,Label='iQScr')

      Call Cho_P_ReoQual(iQScr,IDKVec,nKVec)

      Call mma_deallocate(iQScr)
      Call iCopy(nSym,nKVec,1,nQual,1)

      Call Cho_Timer(C2,W2)
      tDecom(1,4) = tDecom(1,4) + C2 - C1
      tDecom(2,4) = tDecom(2,4) + W2 - W1

!     Compute vectors in each symmetry block.
!     ---------------------------------------

      kV = 1
      kI = 1
      kQD = 1
      Do iSym = 1,nSym

!        Cycle loop if nothing to do in this symmetry.
!        ---------------------------------------------

         If (nQual(iSym) .lt. 1) Go To 100

!        Set vector information.
!        -----------------------

         Call Cho_P_SetVecInf(nQual(iSym),iSym,iPass)

!        Allocate memory for integrals/vectors.
!        --------------------------------------

         l_xInt = max(nnBstR(iSym,2)*nQual(iSym),1)
         Call mma_allocate(xInt,l_xInt,Label='xInt')

         If (nnBstR(iSym,2) .gt. 0) Then

!           Read integral columns from disk, ordered according to IDK.
!           ----------------------------------------------------------

            Call Cho_Timer(C1,W1)
            Call Cho_RdQCol_Indx(xInt,IDKVec(kI),nnBstR(iSym,2),
     &                           nQual(iSym),LuSel(iSym))
            Call Cho_Timer(C2,W2)
            tDecom(1,1) = tDecom(1,1) + C2 - C1
            tDecom(2,1) = tDecom(2,1) + W2 - W1

!           Compute vectors.
!           ----------------

            Call mma_maxDBLE(l_Wrk1)
            Call mma_allocate(Wrk1,l_Wrk1,Label='Wrk1')

            Call Cho_CompVec(Diag,xInt,KVec(kV),QDiag(kQD),
     &                       Wrk1,SIZE(Wrk1),iSym,iPass)

            Call mma_deallocate(Wrk1)

!           Write vectors to disk and update vector counters.
!           -------------------------------------------------

            Call Cho_Timer(C1,W1)
            iVec1 = NumCho(iSym) + 1
            Call Cho_PutVec(xInt,nnBstR(iSym,2),nQual(iSym),
     &                      iVec1,iSym)
            Call Cho_VecBuf_Copy(xInt,nQual(iSym),iSym)
            NumCho(iSym) = NumCho(iSym) + nQual(iSym)
            NumChT = NumChT + nQual(iSym)
            Call Cho_Timer(C2,W2)
            tDecom(1,2) = tDecom(1,2) + C2 - C1
            tDecom(2,2) = tDecom(2,2) + W2 - W1

         End If

!        Transpose vectors on disk (parallel only).
!        ------------------------------------------

         Call Cho_Timer(C1,W1)
         iRed = 2
         Jin = NumV(iSym) + 1
         Jfi = NumV(iSym) + nQual(iSym)
         Call Cho_P_VecTransp(xInt,Jin,Jfi,iSym,iRed,iPass)
         Call Cho_Timer(C2,W2)
         tDecom(1,2) = tDecom(1,2) + C2 - C1
         tDecom(2,2) = tDecom(2,2) + W2 - W1

!        Deallocate memory for integrals/vectors.
!        ----------------------------------------

         Call mma_deallocate(xInt)

!        Empty symmetry blocks jump here.
!        --------------------------------

  100    Continue
         kV = kV + nQual(iSym)**2
         kI = kI + nQual_Old(iSym)
         kQD = kQD + nQual_Old(iSym)

      End Do

!     Deallocations.
!     --------------

      Call mma_deallocate(LQ_Tot)
      Call mma_deallocate(QDiag)
      Call mma_deallocate(IDKVec)
      Call mma_deallocate(KVec)

!     Print.
!     ------

      If (iPrint .ge. Inf_Progress) Then
         Do iSym = 1,nSym
            NumCho_Old(iSym) = NumCho(iSym) - NumCho_Old(iSym)
         End Do
         Write(Lupri,'(80A)') ('-',I=1,LenLin)
         Write(Lupri,'(A,8I8)')
     &   '#vec. gener.  : ',(NumCho_OLD(iSym),iSym=1,nSym)
      Else If (iPrint .ge. Inf_Pass) Then
         Do iSym = 1,nSym
            NumCho_Old(iSym) = NumCho(iSym) - NumCho_Old(iSym)
         End Do
         Write(Lupri,'(A,8I8)')
     &   '#vec. gener.  : ',(NumCho_OLD(iSym),iSym=1,nSym)
      End If

      End
