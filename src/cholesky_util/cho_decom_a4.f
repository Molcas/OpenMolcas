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
      Implicit Real*8 (a-h,o-z)
      Real*8  Diag(*)
      Integer LstQSP(NumSP)
#include "cholesky.fh"
#include "chovecbuf.fh"
#include "choptr.fh"
#include "cholq.fh"
#include "choprint.fh"
#include "WrkSpc.fh"

      Character*12 SecNam
      Parameter (SecNam = 'Cho_Decom_A4')

      Integer NumCho_Old(8), nQual_Old(8)
      Integer NumV(8)

      nKVec(i)=iWork(ip_nKVec-1+i)

#if defined (_DEBUG_)
#endif

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
      l_LQ_Sym(1) = nQual(1)*NumV(1)
      l_LQ = l_LQ_Sym(1)
      Do iSym = 2,nSym
         l_KVec = l_KVec + nQual(iSym)**2
         l_IDKVec = l_IDKVec + nQual(iSym)
         l_LQ_Sym(iSym) = nQual(iSym)*NumV(iSym)
         l_LQ = l_LQ + l_LQ_Sym(iSym)
      End Do
      l_QDiag = l_IDKVec
      l_lQ = max(l_LQ,1) ! because there might not be any prev. vecs.
      Call GetMem('KVec','Allo','Real',ip_KVec,l_KVec)
      Call GetMem('IDKVec','Allo','Inte',ip_IDKVec,l_IDKVec)
      Call GetMem('QDiag','Allo','Real',ip_QDiag,l_QDiag)

C     Extract elements corresponding to qualified diagonals from
C     previous Cholesky vectors (if any).
C     ----------------------------------------------------------

      Call Cho_Timer(C1,W1)
      Call GetMem('LQ','Allo','Real',ip_LQ,l_LQ)
      ip_LQ_Sym(1) = ip_LQ
      Do iSym = 2,nSym
         ip_LQ_Sym(iSym) = ip_LQ_Sym(iSym-1) + l_LQ_Sym(iSym-1)
      End Do
      Call Cho_P_GetLQ(Work(ip_LQ),l_LQ,LstQSP,NumSP)
      Call Cho_Timer(C2,W2)
      tDecom(1,2) = tDecom(1,2) + C2 - C1
      tDecom(2,2) = tDecom(2,2) + W2 - W1

C     Extract qualified diagonal integral block.
C     ------------------------------------------

      Call Cho_Timer(C1,W1)

      l_nKVec = nSym
      l_MQ = l_KVec
      Call GetMem('nKVec','Allo','Inte',ip_nKVec,l_nKVec)
      Call GetMem('MQ','Allo','Real',ip_MQ,l_MQ)
      Call Cho_P_GetMQ(Work(ip_MQ),l_MQ,LstQSP,NumSP)

C     Decompose qualified diagonal block.
C     The qualified diagonals are returned in QDiag.
C     ----------------------------------------------

      Call Cho_Dec_Qual(Diag,Work(ip_LQ),Work(ip_MQ),Work(ip_Kvec),
     &                  iWork(ip_IDKVec),iWork(ip_nKVec),
     &                  Work(ip_QDiag))

C     Deallocate MQ.
C     --------------

      Call GetMem('MQ','Free','Real',ip_MQ,l_MQ)

C     Reorder the elements of the K-vectors according to IDK ordering.
C     ----------------------------------------------------------------

      MxQ = nQual(1)
      Do iSym = 2,nSym
         MxQ = max(MxQ,nQual(iSym))
      End Do

      l_KVScr = MxQ
      Call GetMem('KVScr','Allo','Real',ip_KVScr,l_KVScr)

      kK1 = ip_KVec - 1
      kK2 = kK1
      kID = ip_IDKVec - 1
      kS  = ip_KVScr - 1
      Do iSym = 1,nSym
         Do iK = 1,nKVec(iSym)
            kK_1 = kK1 + nQual(iSym)*(iK-1) + 1
            Call dCopy_(nQual(iSym),Work(kK_1),1,Work(ip_KVScr),1)
            kK_2 = kK2 + nKVec(iSym)*(iK-1)
            Do jK = 1,nKVec(iSym)
               lK = iWork(kID+jK)
               Work(kK_2+jK) = Work(kS+lK)
            End Do
         End Do
         kK1 = kK1 + nQual(iSym)**2
         kK2 = kK2 + nKVec(iSym)**2
         kID = kID + nQual(iSym)
      End Do

C     Reorder QDiag to IDK ordering.
C     ------------------------------

      kID = ip_IDKVec - 1
      kQD = ip_QDiag - 1
      Do iSym = 1,nSym
         Call dCopy_(nQual(iSym),Work(kQD+1),1,Work(ip_KVScr),1)
         Do iK = 1,nKVec(iSym)
            lK = iWork(kID+iK)
            Work(kQD+iK) = Work(ip_KVScr-1+lK)
         End Do
         kQD = kQD + nQual(iSym)
         kID = kID + nQual(iSym)
      End Do

C     Reorder elements of LQ vectors to IDK ordering.
C     -----------------------------------------------

      kID = ip_IDKVec - 1
      Do iSym = 1,nSym
         Do jVec = 1,NumV(iSym)
            kLQ = ip_LQ_Sym(iSym) + nQual(iSym)*(jVec-1) - 1
            Call dCopy_(nQual(iSym),Work(kLQ+1),1,Work(ip_KVScr),1)
            Do iK = 1,nKVec(iSym)
               lK = iWork(kID+iK)
               Work(kLQ+iK) = Work(ip_KVScr-1+lK)
            End Do
         End Do
         kID = kID + nQual(iSym)
         ldLQ(iSym) = max(nQual(iSym),1) ! leading dim of LQ
      End Do

      Call GetMem('KVScr','Free','Real',ip_KVScr,l_KVScr)

C     Reset qualification index arrays to IDK ordering.
C     Local as well as global are reordered.
C     -------------------------------------------------

      Call iCopy(nSym,nQual,1,nQual_Old,1)
      l_iQScr = MxQ
      Call GetMem('iQScr','Allo','Inte',ip_iQScr,l_iQScr)
      Call Cho_P_ReoQual(iWork(ip_iQScr),iWork(ip_IDKVec),
     &                   iWork(ip_nKVec))
      Call GetMem('iQScr','Free','Inte',ip_iQScr,l_iQScr)
      Call iCopy(nSym,iWork(ip_nKVec),1,nQual,1)

C     Deallocate K-vector counter (now stored as nQual counter).
C     ----------------------------------------------------------

      Call GetMem('nKVec','Free','Inte',ip_nKVec,l_nKVec)

      Call Cho_Timer(C2,W2)
      tDecom(1,4) = tDecom(1,4) + C2 - C1
      tDecom(2,4) = tDecom(2,4) + W2 - W1

C     Compute vectors in each symmetry block.
C     ---------------------------------------

      kV = ip_KVec
      kI = ip_IDKVec
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
            Call Cho_RdQCol_Indx(Work(ip_xInt),iWork(kI),
     &                           nnBstR(iSym,2),nQual(iSym),
     &                           LuSel(iSym))
            Call Cho_Timer(C2,W2)
            tDecom(1,1) = tDecom(1,1) + C2 - C1
            tDecom(2,1) = tDecom(2,1) + W2 - W1

C           Compute vectors.
C           ----------------

            Call GetMem('CmpV_Max','Max ','Real',ip_Wrk1,l_Wrk1)
            Call GetMem('CmpV_Wrk','Allo','Real',ip_Wrk1,l_Wrk1)
            Call Cho_CompVec(Diag,Work(ip_xInt),Work(kV),Work(kQD),
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

      Call GetMem('LQ','Free','Real',ip_LQ,l_LQ)
      ip_LQ = -999999
      l_LQ = 0
      Do iSym = 1,nSym
         ip_LQ_Sym(iSym) = -999999
         l_LQ_Sym(iSym) = 0
         ldLQ(iSym) = 0
      End Do
      Call GetMem('QDiag','Free','Real',ip_QDiag,l_QDiag)
      Call GetMem('IDKVec','Free','Inte',ip_IDKVec,l_IDKVec)
      Call GetMem('KVec','Free','Real',ip_KVec,l_KVec)

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

#if defined (_DEBUG_)
#endif

      End
