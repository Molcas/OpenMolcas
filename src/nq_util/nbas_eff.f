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
      Integer Function nBas_Eff(NrExp,NrBas,Exp,Cff,nExp_Eff)
      Implicit Real*8 (a-h,o-z)
      Real*8 Exp(NrExp), Cff(NrExp,NrBas)
*
      nBas_Eff=NrBas
*
      Do iBas = 1, NrBas
*
         Do iExp = 1, nExp_Eff
*
            If (Cff(iExp,iBas).ne.0.0D0) Then
               nBas_Eff = NrBas-iBas+1
               Return
            End If
*
         End Do
*
      End Do
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(Exp)
      End
