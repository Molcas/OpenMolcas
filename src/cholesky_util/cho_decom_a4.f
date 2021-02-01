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
      SubRoutine Cho_Decom_A4(Diag,LstQSP,NumSP,iPass)
C
C     Purpose: decompose qualified columns ("parallel" algorithm).
C
      use ChoArr, only: LQ_Tot, LQ
      use ChoVecBuf, only: nVec_in_Buf

      Implicit Real*8 (a-h,o-z)
      Real*8  Diag(*)
      Integer LstQSP(NumSP)
#include "cholesky.fh"
#include "choprint.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"

      Character(LEN=12), Parameter:: SecNam = 'Cho_Decom_A4'

      Integer NumCho_Old(8), nQual_Old(8)
      Integer NumV(8), nkVec(8)

      Real*8, Allocatable:: KVScr(:), MQ(:), KVec(:)
      Integer, Allocatable:: IDKVec(:)
*                                                                      *
************************************************************************
*                                                                      *
      Interface
      SubRoutine Cho_P_GetLQ(QVec,l_QVec,LstQSP,nQSP)
      Integer l_QVec, nQSP
      Real*8, Target::  QVec(l_Qvec)
      Integer LstQSP(nQSP)
      End SubRoutine Cho_P_GetLQ
      End Interface
*                                                                      *
************************************************************************
*                                                                      *

C     Print header.
C     -------------

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

C     Allocations.
C     ------------

      Call Cho_P_GetGV(numV,nSym)

      l_KVec = nQual(1)**2
      l_IDKVec = nQual(1)
      l_LQ = nQual(1)*NumV(1)
      Do iSym = 2,nSym
         l_KVec = l_KVec + nQual(iSym)**2
         l_IDKVec = l_IDKVec + nQual(iSym)
         l_LQ = l_LQ + nQual(iSym)*NumV(iSym)
      End Do
      l_QDiag = l_IDKVec
      l_LQ = max(l_LQ,1) ! because there might not be any prev. vecs.
      Call mma_allocate(KVec,l_KVec,Label='KVec')
      Call mma_allocate(IDKVec,l_IDKVec,Label='IDKVec')
      Call GetMem('QDiag','Allo','Real',ip_QDiag,l_QDiag)

C     Extract elements corresponding to qualified diagonals from
C     previous Cholesky vectors (if any).
C     ----------------------------------------------------------

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

C     Extract qualified diagonal integral block.
C     ------------------------------------------

      Call Cho_Timer(C1,W1)

      Call mma_allocate(MQ,l_KVec,Label='MQ')
      Call Cho_P_GetMQ(MQ,SIZE(MQ),LstQSP,NumSP)

C     Decompose qualified diagonal block.
C     The qualified diagonals are returned in QDiag.
C     ----------------------------------------------

      Call Cho_Dec_Qual(Diag,LQ_Tot,MQ,KVec,
     &                  IDKVec,nKVec,Work(ip_QDiag))

C     Deallocate MQ.
C     --------------

      Call mma_deallocate(MQ)

C     Reorder the elements of the K-vectors according to IDK ordering.
C     ----------------------------------------------------------------

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

C     Reorder QDiag to IDK ordering.
C     ------------------------------

      kID = 0
      kQD = ip_QDiag - 1
      Do iSym = 1,nSym
         Call dCopy_(nQual(iSym),Work(kQD+1),1,KVScr,1)
         Do iK = 1,nKVec(iSym)
            lK = IDKVec(kID+iK)
            Work(kQD+iK) = KVScr(lK)
         End Do
         kQD = kQD + nQual(iSym)
         kID = kID + nQual(iSym)
      End Do

C     Reorder elements of LQ vectors to IDK ordering.
C     -----------------------------------------------

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

C     Reset qualification index arrays to IDK ordering.
C     Local as well as global are reordered.
C     -------------------------------------------------

      Call iCopy(nSym,nQual,1,nQual_Old,1)
      l_iQScr = MxQ
      Call GetMem('iQScr','Allo','Inte',ip_iQScr,l_iQScr)
      Call Cho_P_ReoQual(iWork(ip_iQScr),IDKVec,nKVec)
      Call GetMem('iQScr','Free','Inte',ip_iQScr,l_iQScr)
      Call iCopy(nSym,nKVec,1,nQual,1)

      Call Cho_Timer(C2,W2)
      tDecom(1,4) = tDecom(1,4) + C2 - C1
      tDecom(2,4) = tDecom(2,4) + W2 - W1

C     Compute vectors in each symmetry block.
C     ---------------------------------------

      kV = 1
      kI = 1
      kQD = ip_QDiag
      Do iSym = 1,nSym

C        Cycle loop if nothing to do in this symmetry.
C        ---------------------------------------------

         If (nQual(iSym) .lt. 1) Go To 100

C        Set vector information.
C        -----------------------

         Call Cho_P_SetVecInf(nQual(iSym),iSym,iPass)

C        Allocate memory for integrals/vectors.
C        --------------------------------------

         l_xInt = max(nnBstR(iSym,2)*nQual(iSym),1)
         Call GetMem('xInt','Allo','Real',ip_xInt,l_xInt)

         If (nnBstR(iSym,2) .gt. 0) Then

C           Read integral columns from disk, ordered according to IDK.
C           ----------------------------------------------------------

            Call Cho_Timer(C1,W1)
            Call Cho_RdQCol_Indx(Work(ip_xInt),IDKVec(kI),
     &                           nnBstR(iSym,2),nQual(iSym),
     &                           LuSel(iSym))
            Call Cho_Timer(C2,W2)
            tDecom(1,1) = tDecom(1,1) + C2 - C1
            tDecom(2,1) = tDecom(2,1) + W2 - W1

C           Compute vectors.
C           ----------------

            Call GetMem('CmpV_Max','Max ','Real',ip_Wrk1,l_Wrk1)
            Call GetMem('CmpV_Wrk','Allo','Real',ip_Wrk1,l_Wrk1)
            Call Cho_CompVec(Diag,Work(ip_xInt),KVec(kV),Work(kQD),
     &                       Work(ip_Wrk1),l_Wrk1,iSym,iPass)
            Call GetMem('CmpV_Wrk','Free','Real',ip_Wrk1,l_Wrk1)

C           Write vectors to disk and update vector counters.
C           -------------------------------------------------

            Call Cho_Timer(C1,W1)
            iVec1 = NumCho(iSym) + 1
            Call Cho_PutVec(Work(ip_xInt),nnBstR(iSym,2),nQual(iSym),
     &                      iVec1,iSym)
            Call Cho_VecBuf_Copy(Work(ip_xInt),nQual(iSym),iSym)
            NumCho(iSym) = NumCho(iSym) + nQual(iSym)
            NumChT = NumChT + nQual(iSym)
            Call Cho_Timer(C2,W2)
            tDecom(1,2) = tDecom(1,2) + C2 - C1
            tDecom(2,2) = tDecom(2,2) + W2 - W1

         End If

C        Transpose vectors on disk (parallel only).
C        ------------------------------------------

         Call Cho_Timer(C1,W1)
         iRed = 2
         Jin = NumV(iSym) + 1
         Jfi = NumV(iSym) + nQual(iSym)
         Call Cho_P_VecTransp(Work(ip_xInt),Jin,Jfi,iSym,iRed,iPass)
         Call Cho_Timer(C2,W2)
         tDecom(1,2) = tDecom(1,2) + C2 - C1
         tDecom(2,2) = tDecom(2,2) + W2 - W1

C        Deallocate memory for integrals/vectors.
C        ----------------------------------------

         Call GetMem('xInt','Free','Real',ip_xInt,l_xInt)

C        Empty symmetry blocks jump here.
C        --------------------------------

  100    Continue
         kV = kV + nQual(iSym)**2
         kI = kI + nQual_Old(iSym)
         kQD = kQD + nQual_Old(iSym)

      End Do

C     Deallocations.
C     --------------

      Call mma_deallocate(LQ_Tot)
      Call GetMem('QDiag','Free','Real',ip_QDiag,l_QDiag)
      Call mma_deallocate(IDKVec)
      Call mma_deallocate(KVec)

C     Print.
C     ------

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
