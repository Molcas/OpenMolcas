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
      Logical Function R_Stab_A(R,S,nS)
      Integer S(nS), R
*
      R_Stab_A = .False.
      Do iS = 1, nS
         If (R.eq.S(iS)) Then
            R_Stab_A = .True.
            Return
         End If
      End Do
*
      Return
      End
