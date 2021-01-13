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
      use Basis_Info, only: dbsc
#include "Basis_Mode_Parameters.fh"
#include "Basis_Mode.fh"
*
      If (dbsc(i)%Aux) Then
         Basis_Mode = Auxiliary_Mode
      Else
         Basis_Mode = Valence_Mode
      End If
*
      Do k = i+1, j
         If (dbsc(i)%Aux.neqv.dbsc(k)%Aux) Then
            Call WarningMessage(2,
     &              'dbsc(i)%Aux.ne.dbsc(k)%Aux')
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
