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
      SubRoutine Cho_P_GetQD(QD)
C
C     Purpose: copy qualified diagonal elements from global diagonal to
C              array QD.
C
      use ChoSwp, only: iQuAB, IndRed_G
      Implicit None
      Real*8 QD(*)
#include "cholesky.fh"
#include "choptr.fh"
#include "choglob.fh"
#include "WrkSpc.fh"

      Integer kQD, kD, iSym, iQ, iAB

      kQD = 0
      kD = ip_Diag_G - 1
      Do iSym = 1,nSym
         Do iQ = 1,nQual(iSym)
            iAB = IndRed_G(iQuAB(iQ,iSym),2)
            QD(kQD+iQ) = Work(kD+iAB)
         End Do
         kQD = kQD + nQual(iSym)
      End Do

      End
