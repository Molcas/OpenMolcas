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
      SubRoutine Cho_VecBuf_Subtr(xInt,Wrk,lWrk,iSym,DoTime,DoStat)
C
C     Purpose: subtract contributions to qualified columns from the
C              vectors stored in the buffer (if any).
C
C     DoTime: time as vector subtraction.
C     DpStat: update statistics info (#calls to dGeMM).
C
      use ChoSwp, only: iQuAB, nnBstRSh, iiBstRSh
      use ChoArr, only: LQ
      use ChoVecBuf
      use ChoSubScr, only: Cho_SScreen, SSTau, SubScrStat, DSubScr,
     &                     ip_DSPNm, SSNorm
#include "implicit.fh"
      Real*8, Target::  xInt(*), Wrk(lWrk)
      Logical DoTime, DoStat
#include "cholesky.fh"
#include "WrkSpc.fh"

      Character(LEN=16), Parameter:: SecNam = 'Cho_VecBuf_Subtr'

#if defined (_DEBUGPRINT_)
      Logical, Parameter:: LocDbg = .true.
#else
      Logical, Parameter:: LocDbg = .false.
#endif

      Real*8, Parameter:: xMOne = -1.0d0, One = 1.0d0

      Real*8, Pointer:: V(:,:)=>Null(), U(:,:)=>Null(), W(:,:)=>Null()

      DSPNm(i)=Work(ip_DSPNm-1+i)

C     Return if nothing to do.
C     ------------------------

      If (l_ChVBuf_Sym(iSym) .lt. 1) Then
         If (LocDbg) Then
            Write(Lupri,*) SecNam,': returns immediately!'
            Write(Lupri,*) ' -- no buffer allocated for sym. ',iSym
         End If
         Return
      End If
      If (nVec_in_Buf(iSym) .lt. 1) Then
         If (LocDbg) Then
            Write(Lupri,*) SecNam,': returns immediately!'
            Write(Lupri,*) ' -- buffer is empty for sym. ',iSym
         End If
         Return
      End If
      If (nQual(iSym) .lt. 1) Then
         If (LocDbg) Then
            Write(Lupri,*) SecNam,': returns immediately!'
            Write(Lupri,*) ' -- no qualified columns of sym. ',iSym
         End If
         Return
      End If
      If (nnBstR(iSym,2) .lt. 1) Then
         If (LocDbg) Then
            Write(Lupri,*) SecNam,': returns immediately!'
            Write(Lupri,*) ' -- empty symmetry block (sym. ',iSym,')'
         End If
         Return
      End If

C     Start timing.
C     -------------

      If (DoTime) Call Cho_Timer(C1,W1)

C     Initialize.
C     -----------

      xTot = 0.0d0
      xDon = 0.0d0

C     Set up vector batch.
C     --------------------

      nVec = min(lWrk/nQual(iSym),nVec_in_Buf(iSym))
      If (nVec .lt. 1) Then
         Call Cho_Quit('Insufficient memory for batch in '//SecNam,101)
         nBatch = -999999 ! avoid compiler warnings
      Else
         nBatch = (nVec_in_Buf(iSym)-1)/nVec + 1
      End If

C     Map the integral array, xInt, onto the pointer U
C     ------------------------------------------------

      lRow = nnBstR(iSym,2)
      lCol = nQual(iSym)
      iS = 1
      iE = iS - 1 + lRow*lCol

      U(1:lRow,1:lCol) => xInt(iS:iE)

C     Start batch loop.
C     -----------------

      Do iBatch = 1,nBatch

C        Set info for this batch.
C        ------------------------

         If (iBatch .eq. nBatch) Then
            NumV = nVec_in_Buf(iSym) - nVec*(nBatch-1)
         Else
            NumV = nVec
         End If
         iVec0 = nVec*(iBatch-1)

#if defined (_DEBUGPRINT_)
         Need = nQual(iSym)*NumV
         If (lWrk .lt. Need) Then
            Call Cho_Quit('Batch setup error in '//SecNam,104)
         End If
#endif

         lRow = nnBstR(iSym,2)
         lCol = iVec0 + NumV
         iS   = ip_ChVBuf_Sym(iSym)
         iE   = iS - 1 + lRow*lCol

         V(1:lRow,1:lCol) => CHVBUF(iS:iE)

C        Screened or unscreened subtraction section.
C        The screened version uses level 2 blas, while the unscreened
C        one employs level 3 blas.
C        ------------------------------------------------------------

         If (Cho_SScreen) Then

            lRow = NumV
            lCol = nQual(iSym)
            iS   = 1
            iE   = iS - 1 + lRow*lCol

            W(1:lRow,1:lCol) => Wrk(iS:iE)

C           Copy out sub-blocks corresponding to qualified diagonals:
C           L(#J,{ab})
C           ---------------------------------------------------------

            Do jVec = 1,NumV
               Do iAB = 1,nQual(iSym)
                  jAB = iQuAB(iAB,iSym) - iiBstR(iSym,2)
                  W(jVec,iAB) = V(jAB,jVec+iVec0)
               End Do
            End Do

C           Subtract:
C           (gd|{ab}) <- (gd|{ab}) - sum_J L(gd,#J) * L(#J,{ab})
C           for each ab in {ab}.
C           ----------------------------------------------------

            Call Cho_SubScr_Dia(V(1,iVec0+1),NumV,iSym,2,SSNorm)

            Do iAB = 1,nQual(iSym)
               Do iShGD = 1,nnShl
                  nGD = nnBstRSh(iSym,iShGD,2)
                  If (nGD < 1 ) Cycle
                  iGD = iiBstRSh(iSym,iShGD,2)
                  xTot = xTot + 1.0d0
                  jAB = iQuab(iAB,iSym) - iiBstR(iSym,2)
                  Tst = sqrt(DSPNm(iShGD)*DSubScr(jAB))
                  If (Tst<=SSTau) Cycle
                  xDon = xDon + 1.0d0
                  Call dGeMV_('N',nGD,NumV,
     &                       xMOne,V(1+iGD,iVec0+1),nnBstR(iSym,2),
     &                       W(1,iAB),1,One,U(1+iGD,iAB),1)
               End Do
            End Do

         Else ! unscreened subtraction

            If (Associated(LQ(iSym)%Array)) Then

C              If the qualified block, L({ab},#J), is already in core,
C              use this block.
C              -------------------------------------------------------


               Call DGEMM_('N','T',nnBstR(iSym,2),nQual(iSym),NumV,
     &                    xMOne,V(1,iVec0+1),nnBstR(iSym,2),
     &                          LQ(iSym)%Array(:,iVec0+1),
     &                          SIZE(LQ(iSym)%Array,1),
     &                    One,U,nnBstR(iSym,2))

            Else

               lRow = nQual(iSym)
               lCol = NumV
               iS   = 1
               iE   = iS - 1 + lRow*lCol

               W(1:lRow,1:lCol) => WrK(iS:iE)

C              Copy out sub-blocks corresponding to qualified diagonals:
C              L({ab},#J).
C              ---------------------------------------------------------

               Do jVec = 1,NumV
                  Do iAB = 1,nQual(iSym)
                     jAB = iQuAB(iAB,iSym) - iiBstR(iSym,2)
                     W(iAB,jVec) = V(jAB,jVec+iVec0)
                  End Do
               End Do

C              Subtract:
C              (gd|{ab}) <- (gd|{ab}) - sum_J L(gd,#J) * L({ab},#J)
C              ----------------------------------------------------


               Call DGEMM_('N','T',nnBstR(iSym,2),nQual(iSym),NumV,
     &                    xMOne,V(1,iVec0+1),nnBstR(iSym,2),
     &                          W,nQual(iSym),
     &                    One,U,nnBstR(iSym,2))

            End If

         End If

         V=>Null()
         U=>Null()
         W=>Null()

      End Do

C     Update statistics info.
C     -----------------------

      If (DoStat) nDGM_Call = nDGM_Call + nBatch
      If (Cho_SScreen) Then
         SubScrStat(1) = SubScrStat(1) + xTot
         SubScrStat(2) = SubScrStat(2) + xDon
      End If

C     Update global timing.
C     ---------------------

      If (DoTime) Then
         Call Cho_Timer(C2,W2)
         tDecom(1,3) = tDecom(1,3) + C2 - C1
         tDecom(2,3) = tDecom(2,3) + W2 - W1
      End If

      End
