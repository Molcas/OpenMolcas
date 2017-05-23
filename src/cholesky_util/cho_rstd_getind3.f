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
      SubRoutine Cho_RstD_GetInd3(iSP2F,l_iSP2F)
      Implicit None
      Integer l_iSP2F
      Integer iSP2F(l_iSP2F)
#include "cholesky.fh"

      Integer iOpt, iAdr

      iOpt = 2
      iAdr = nSym*nnShl + 2*nnBstRT(1)
      Call iDAFile(LuRed,iOpt,iSP2F,l_iSP2F,iAdr)

      End
