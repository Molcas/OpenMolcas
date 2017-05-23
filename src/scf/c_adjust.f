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
      Subroutine  C_Adjust(CInter,n,CThr)
      Implicit None
      Integer n, i
      Real*8 CInter(n), CThr, Fact
*
      If (CInter(n).lt.CThr) Then
         Fact=(1.0D0-CThr)/(1.0D0-CInter(n))
         Do i = 1, n - 1
            CInter(i) = Fact * CInter(i)
         End Do
         CInter(n)=CThr
      End If
*
      Return
      End
