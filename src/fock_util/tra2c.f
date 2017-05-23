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
      Subroutine Tra2c(i,iSym,iBas,iAsh,
     &                 j,jSym,jBas,jAsh,
     &                 kl_Bas_pairs,ij_Orb_pairs,
     &                 CMO_i,CMO_j,IJKL,C,TUKL)
      Real*8 CMO_i(iBas,iAsh), CMO_j(jBas,jAsh),
     &       C(ij_Orb_pairs), IJKL(kl_Bas_pairs),
     &       TUKL(kl_Bas_pairs,ij_Orb_pairs)
*
*     Call RecPrt('IJKL',' ',IJKL,1,kl_Bas_Pairs)
*
      If (iSym.eq.jSym.and.i.ne.j) Then
         ijA = 0
         Do iA = 1, iAsh
            Do jA = 1, iA
               ijA = ijA + 1
               C(ijA)=CMO_i(i,iA)*CMO_i(j,jA)
     &               +CMO_i(j,iA)*CMO_i(i,jA)
            End Do
         End Do
*        Call TriPrt('C',' ',C,iAsh)
      ElseIf (iSym.eq.jSym) Then
         ijA = 0
         Do iA = 1, iAsh
            Do jA = 1, iA
               ijA = ijA + 1
               C(ijA)=CMO_i(i,iA)*CMO_i(i,jA)
            End Do
         End Do
*        Call TriPrt('C',' ',C,iAsh)
      Else
         ijA = 0
         Do iA = 1, iAsh
            Do jA = 1, jAsh
               ijA = ijA + 1
               C(ijA)=CMO_i(i,iA)*CMO_j(j,jA)
            End Do
         End Do
*        Call RecPrt('C',' ',C,jAsh,iAsh)
      End If
*     Write(6,*) ' ij_Orb_Pairs=',ij_Orb_Pairs
*
      Call DNaXpY(ij_Orb_pairs,kl_Bas_pairs,
     &            C,1,
     &            IJKL,1,0,
     &            TUKL,1,kl_Bas_pairs)
*
      Return
      End
