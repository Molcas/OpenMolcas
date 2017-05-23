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
      Subroutine  ReMap_U_k(U_k,nU_k,U_k_New,nU_k_New,iSO_ab)
      Implicit Real*8 (A-H,O-Z)
      Real*8 U_k(nU_k), U_k_New(nU_k_New)
      Integer iSO_ab(2,nU_k)
*
      Do k=1,nU_k
         i=iSO_ab(1,k)
         j=iSO_ab(2,k)
         ij=i*(i-1)/2 + j
         If (i.eq.j) Then
            U_k_New(ij) = U_k(k)
         Else
            U_k_New(ij) = 0.5D0*U_k(k)
         End If
      End Do

      Return
      End
