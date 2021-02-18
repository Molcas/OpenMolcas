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
      SUBROUTINE Cho_P_QualSwp()
      use ChoSwp, only: iQuAB, iQuAB_L, pTemp
      use ChoArr, only: nQual_L
      Implicit None
#include "cholesky.fh"

      Integer i, scr

C --- Swap nQual array  local <-> global
      Do i=1,nSym
         scr = nQual_L(i)
         nQual_L(i) = nQual(i)
         nQual(i) = scr
      End Do

C --- Swap pointer for iQuAB  local <-> global
      pTemp => iQuAB
      iQuAB => iQuAB_L
      iQuAB_L => pTemp

      End
