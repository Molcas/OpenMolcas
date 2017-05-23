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
      SubRoutine Cho_ReoQual(iQuAB,MxQ,nSym,iQScr,IDK,nK,nQ)
C
C     Purpose: reorder iQuAB array to IDK ordering.
C
      Implicit None
      Integer MxQ, nSym
      Integer iQuAB(MxQ,nSym)
      Integer iQScr(MxQ)
      Integer IDK(*)
      Integer nK(nSym)
      Integer nQ(nSym)

      Integer iSym, kID
      Integer iK, lK

      kID = 0
      Do iSym = 1,nSym
         If (nQ(iSym) .gt. 0) Then
            Call iCopy(nQ(iSym),iQuAB(1,iSym),1,iQScr,1)
            Do iK = 1,nK(iSym)
               lK = IDK(kID+iK)
               iQuAB(iK,iSym) = iQScr(lK)
            End Do
            kID = kID + nQ(iSym)
         Else
            Do iK = 1,nK(iSym)
               iQuAB(iK,iSym) = 0
            End Do
         End If
      End Do

      End
