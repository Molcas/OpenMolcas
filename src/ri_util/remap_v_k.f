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
      Subroutine  ReMap_V_k(iSym,V_k,nV_k,V_k_New,nV_k_New,iSO_ab,ij2K,
     &                      m_ij2K)
      Implicit Real*8 (A-H,O-Z)
      Real*8 V_k(nV_k), V_k_New(nV_k_New)
      Integer iSym, iSO_ab(2,nV_k), ij2K(m_ij2K)
*
      If (iSym .eq. 1) Then
         Do k=1,nV_k
            i=iSO_ab(1,k)
            j=iSO_ab(2,k)
            ij=i*(i-1)/2 + j
            If (i.eq.j) Then
               V_k_New(ij) = V_k(k)
            Else
               V_k_New(ij) = 0.5D0*V_k(k)
            End If
            ij2K(ij)=k
         End Do
*
c        write(6,*) 'Triang <Vk|Vk> : ',ddot_(nV_k_New,V_k_New,1,
c     &                                               V_k_New,1)
*
      Else

         Do k=1,nV_k
            i=iSO_ab(1,k)
            j=iSO_ab(2,k)
            ij=i*(i-1)/2 + j
            ij2K(ij)=k
         End Do
      EndIf
*
      Return
      End
