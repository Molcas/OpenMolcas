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
      Subroutine BPut(EVec,nDim,BMx,nX,Smmtrc,nQQ,Degen)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 EVec(nDim,nQQ), BMx(nX,nQQ), Degen(nX)
      Logical Smmtrc(nX)
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUG_
*                                                                      *
************************************************************************
*                                                                      *
      iDim=0
      Do iX = 1, nX
         If (Smmtrc(iX)) Then
            iDim=iDim+1
            Do iQQ= 1, nQQ
               BMx(iX,iQQ)=EVec(iDim,iQQ)/Sqrt(Degen(iX))
            End Do
         Else
            Do iQQ = 1, nDim
               BMx(iX,iQQ)=Zero
            End Do
         End If
      End Do
#ifdef _DEBUG_
      Call RecPrt('BPut: BMx',' ',BMx,nX,nQQ)
#endif
*
      Return
      End
