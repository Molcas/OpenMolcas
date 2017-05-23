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
      Subroutine Fix_Coeff(nPrim,nCntrc,Coeff_c,Coeff_p,
     &                     Mode)
      Implicit Real*8 (a-h,o-z)
      Real*8 Coeff_c(nPrim,nCntrc), Coeff_p(nPrim,nPrim)
      Character Mode*1
*
*     Put in the normalization constant for the product
*     basis function.
*
      If (Mode.eq.'F') Then
         Do iC = 1, nCntrc
            Do iP = 1, nPrim
               Coeff_c(iP,iC) = Coeff_c(iP,iC)/Coeff_p(iP,iP)
            End Do
         End Do
      Else
         Do iC = 1, nCntrc
            Do iP = 1, nPrim
               Coeff_c(iP,iC) = Coeff_c(iP,iC)*Coeff_p(iP,iP)
            End Do
         End Do
      End If
*
      Return
      End
