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
      Integer Function iUR(iR,iU)
      iUR=0
      Do i = 0, 7
         If (iAnd(iU,2**i).eq.2**i) Then
            iUR = iOr(iUR,2**iEor(i,iR))
         End If
      End Do
*     Write (*,*) ' iUR=',iUR
      Return
      End
