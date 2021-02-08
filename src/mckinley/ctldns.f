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
        SubRoutine CtlDns(iDCRR,iDCRS,iDCRT,jOp)

      Integer jop(6)

*              Djl. Some care has to be taken here. Assume that there
*              are two operators, T and S which generates the center
*              pairs A,T(B) and A,S(B). If these pairs are symmetry
*              related we will only
*
*--------------Dij
               iR = iDCRR
               jOp(1) = NrOpr(iR) + 1
*--------------Dkl
               iS = iDCRS
               jOp(2) = NrOpr(iS) + 1
*--------------Dik
               iT  = iDCRT
               jOp(3) = NrOpr(iT) + 1
*--------------Dil
               iTS = iEor(iT,iS)
               jOp(4) = NrOpr(iTS) + 1
*--------------Djk
               iRT = iEor(iR,iT)
               jOp(5) = NrOpr(iRT) + 1
*--------------Djl
               iRTS= iEor(iRT,iS)
               jOp(6) = NrOpr(iRTS) + 1
*
      Return
      End
*
