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
      Integer Function iRnge(Val,Bin,nBin)
      Implicit None
      Integer nBin
      Real*8 Val, Bin(nBin)

      Integer iBin

      iRnge = nBin
      Do iBin = 1,nBin-1
         If (Val .gt. Bin(iBin)) Then
            iRnge = iBin
            Return
         End If
      End Do

      End
