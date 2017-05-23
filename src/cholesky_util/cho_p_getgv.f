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
      SUBROUTINE Cho_P_getGV(numV,nSym)

      Implicit None
      Integer nSym, numV(nSym)
#include "cho_para_info.fh"

      If (Cho_Real_Par) Then
         Call Cho_P_getGV_P(numV,nSym)
      Else
         Call Cho_P_getGV_S(numV,nSym)
      End If

      End
*****************************************
      SUBROUTINE Cho_P_getGV_P(numV,nSym)

      Implicit None
      Integer nSym, numV(nSym)

#include "choglob.fh"

      Integer i

      Do i=1,nSym
         numV(i) = numCho_G(i)
      End Do

      End
*****************************************
      SUBROUTINE Cho_P_getGV_S(numV,mSym)

      Implicit None
      Integer mSym, numV(mSym)

#include "cholesky.fh"

      Integer i

      Do i=1,mSym
         numV(i) = numCho(i)
      End Do

      End
