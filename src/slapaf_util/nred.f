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
      Subroutine NRed(ArrIn,ArrOut,nX,nDim,Smmtrc)
      Implicit Real*8 (a-h,o-z)
      Real*8 ArrIn(nX), ArrOut(nDim)
      Logical Smmtrc(nX)
*
      iDim = 0
      Do iX = 1, nX
         If (Smmtrc(iX)) Then
            iDim = iDim + 1
            ArrOut(iDim)=ArrIn(iX)
         End If
      End Do
      If (iDim.ne.nDim) Then
         Write (6,*) 'In NRed: iDim.ne.nDim'
         Call Abend
      End If
*     Call RecPrt('ArrIn',' ',ArrIn,nX,1)
*     Call RecPrt('ArrOut',' ',ArrOut,nDim,1)
*
      Return
      End
