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
      SubRoutine Cho_GetLQ(QVec,l_QVec,LstQSP,nQSP)
C
C     Purpose: extract elements corresponding to qualified columns from
C              the Cholesky vectors in buffer and/or on disk.
C
      use ChoVecBuf, only: nVec_in_Buf
      Implicit None
      Integer l_QVec, nQSP
      Real*8, Target::  QVec(l_Qvec)
      Integer LstQSP(nQSP)
#include "cholesky.fh"

      Character*9 SecNam
      Parameter (SecNam = 'Cho_GetLQ')

      Integer iV1(8), nV(8)
      Integer nTot, iSym

C     Check input.
C     ------------

      If (nQSP .lt. 1) Return
      nTot = NumCho(1)
      Do iSym = 2,nSym
         nTot = nTot + NumCho(iSym)
      End Do
      If (nTot .lt. 1) Return

      nTot = nQual(1)
      Do iSym = 2,nSym
         nTot = nTot + nQual(iSym)
      End Do
      If (nTot .lt. 1) Return

C     Extract from vectors in buffer.
C     -------------------------------

      Call Cho_VecBuf_GetLQ(QVec,l_QVec)

C     Extract from vectors on disk.
C     -----------------------------

      Do iSym = 1,nSym
         iV1(iSym) = nVec_in_Buf(iSym) + 1
         nV(iSym) = NumCho(iSym) - nVec_in_Buf(iSym)
      End Do
      Call Cho_VecDsk_GetLQ(QVec,l_QVec,LstQSP,nQSP,iV1,nV,nSym)

      End
