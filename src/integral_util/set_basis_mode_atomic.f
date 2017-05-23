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
      Subroutine Set_Basis_Mode_Atomic(i,j)
#include "Basis_Mode_Parameters.fh"
#include "Basis_Mode.fh"
#include "itmax.fh"
#include "info.fh"
*
      If (AuxCnttp(i)) Then
         Basis_Mode = Auxiliary_Mode
      Else
         Basis_Mode = Valence_Mode
      End If
*
      Do k = i+1, j
         If (AuxCnttp(i).neqv.AuxCnttp(k)) Then
            Call WarningMessage(2,
     &              'AuxCnttp(i).ne.AuxCnttp(k)')
            Call Abend()
         End If
      End Do
*
      Atomic=.True.
      kCnttp = i
      lCnttp = j
*
      Return
      End
