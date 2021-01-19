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
      SubRoutine Cho_CompVec(Diag,xInt,VecK,QDiag,Wrk,lWrk,iSym,iPass)
C
C     Purpose: compute vectors from decomposition of the qualified
C              diagonal integral block. VecK contains the nQual(iSym)
C              vectors from that decomposition. xInt contains the
C              integral columns. The vectors are returned in xInt.
C              Wrk (dimension lWrk) is work space allocated in calling
C              routine. QDiag contains the qualified diagonals.
C
      use ChoSwp, only: IndRed
      Implicit Real*8 (a-h,o-z)
      Real*8 Diag(*), xInt(*), VecK(*), QDiag(*), Wrk(lWrk)
#include "cholesky.fh"
#include "choptr.fh"
#include "choprint.fh"

      Character*11 SecNam
      Parameter (SecNam = 'Cho_CompVec')

      Logical LocDbg
      Parameter (LocDbg = .False.)

      Integer  Cho_P_IndxParentDiag
      External Cho_P_IndxParentDiag

C     Subtract previous vectors.
C     --------------------------

      Call Cho_Subtr(xInt,Wrk,lWrk,iSym)

C     Debug: check diagonal elements in updated integrals.
C     ----------------------------------------------------

      If (Cho_DiaChk .or. LocDbg) Then
         Tol = Tol_DiaChk
         nErr = 0
         Call Cho_P_ChkInt(xINT,Diag,iSym,nErr,Tol,.True.)
         If (nErr .ne. 0) Then
            Write(Lupri,*) SecNam,': ',nErr,' diagonal errors found!'
            Write(Lupri,*) '          Integral pass: ',iPass
            Write(Lupri,*) '          #Tests: ',nQual(iSym)
c           Write(Lupri,*) '          Printing integrals:'
c           Call Cho_Output(xINT,
c    &                      1,nnBstR(iSym,2),1,nQual(iSym),
c    &                      nnBstR(iSym,2),nQual(iSym),1,Lupri)
            Call Cho_Quit('Diagonal errors in '//SecNam,104)
         Else
            Write(Lupri,*) SecNam,': comparison of qual. integrals ',
     &                  'and current diagonal: no errors !'
         End If
      End If

C     Set max. diagonal for screening.
C     --------------------------------

      QDmax = QDiag(1)

C     Compute vectors.
C     ----------------

      Do i = 1,nQual(iSym)

         kOff0 = nnBstR(iSym,2)*(i-1)
         kK0 = nQual(iSym)*(i-1)

C        Compute vector.
C        ---------------

         xC = QDiag(i)
         Fac = 1.0d0/sqrt(abs(xC))
         kOff = kOff0 + 1
         Call dScal_(nnBstR(iSym,2),Fac,xInt(kOff),1)

C        Zero elements corresponding to zero diagonals.
C        ----------------------------------------------

         Do jAB = 1,nnBstR(iSym,2)
            jAB1 = IndRed(iiBstR(iSym,2)+jAB,2)
            If (Diag(jAB1) .eq. 0.0d0) Then
               xInt(kOff0+jAB) = 0.0d0
            End If
         End Do

C        Update diagonal.
C        ----------------

         Do jAB = 1,nnBstR(iSym,2)
            jAB1 = IndRed(iiBstR(iSym,2)+jAB,2)
            kOff = kOff0 + jAB
            Diag(jAB1) = Diag(jAB1) - xInt(kOff)**2
         End Do

C        Update Qdiag.
C        -------------

         Do j = i,nQual(iSym)
            QDiag(j) = QDiag(j) - VecK(kK0+j)**2
         End Do
         OlDiag = QDiag(i)
         QDiag(i) = 0.0d0

C        Get index (1st reduced set) of the parent diagonal for this
C        vector (global index!).
C        -----------------------------------------------------------

         iABG = Cho_P_IndxParentDiag(i,iSym)

C        Zero treated diagonal element.
C        ------------------------------

         Call Cho_P_ZeroDiag(Diag,iSym,iABG)

C        Check for and zero negative diagonals, count converged, and
C        find max. diagonal element.
C        -----------------------------------------------------------

         Call Cho_ChkDia_A4(Diag,QDmax,iSym,nNeg,nNegT,nConv,xM,yM,zM)

C        nnZTot is the total number of zeroed negative diagonals,
C        updated for statistics purposes.
C        --------------------------------------------------------

         nnZTot = nnZTot + nNeg

C        Subtract this vector from remaining integral columns,
C        M([gd],{ab}_j) -= K({ab}_j,i)*L([gd],NumCho+i)
C        -----------------------------------------------------

         kOff = kOff0 + 1
         Do j = i+1,nQual(iSym)
            Fac = -VecK(kK0+j)
            kInt = nnBstR(iSym,2)*(j-1) + 1
            Call dAXPY_(nnBstR(iSym,2),Fac,xInt(kOff),1,xInt(kInt),1)
         End Do

C        Print.
C        ------

         If (iPrint .ge. Inf_Progress) Then
            iVec = NumCho(iSym) + i
            iVecT = NumChT + i
            Do j = i+1,nQual(iSym)
               xM = max(xM,QDiag(j))
            End Do
            Write(Lupri,'(I3,3(1X,I9),2(1X,D11.3),2(1X,I4),1X,D11.3)')
     &      iSym,iVec,iVecT,iABG,xC,OlDiag,nConv,nNeg,xM
         End If

      End Do

      If (iPrint .ge. Inf_Progress) Call Cho_Flush(Lupri)

      End
