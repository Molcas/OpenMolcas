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

      Implicit None
#include "cholesky.fh"
#include "choptr.fh"
#include "cholq.fh"

      Integer i, scr

C --- Swap nQual array  local <-> global
      Do i=1,nSym
         scr = nQual_L(i)
         nQual_L(i) = nQual(i)
         nQual(i) = scr
      End Do
c      Call iSwap(nSym,nQual,1,nQual_L,1)


C --- Swap pointer for iQuAB  local <-> global
      scr = ip_iQuAB_L
      ip_iQuAB_L = ip_iQuAB
      ip_iQuAB = scr

C --- Swap length for iQuAB   local <-> global
      scr = l_iQuAB_L
      l_iQuAB_L = l_iQuAB
      l_iQuAB = scr


      End
