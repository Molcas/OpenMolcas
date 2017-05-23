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
      SubRoutine Cho_SubScr_Init()
C
C     Purpose: initialize screening in vector subtraction.
C
#include "implicit.fh"
#include "cholesky.fh"
#include "chosubscr.fh"
#include "WrkSpc.fh"

      l_DSubScr = nnBstR(1,1)
      Do iSym = 2,nSym
         l_DSubScr = max(l_DSubScr,nnBstR(iSym,1))
      End Do
      Call Cho_Mem('DSubScr','Allo','Real',ip_DSubScr,l_DSubScr)

      l_DSPNm = nnShl
      Call Cho_Mem('DSPMx','Allo','Real',ip_DSPNm,l_DSPNm)

      End
