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
      Function Multipole_Expansion(q_A,q_B,Dipole_A,Dipole_B,
     &                  R_A,R_B,R,Absolute)
      Implicit Real*8 (A-H,O-Z)
      Real*8 Multipole_Expansion
      Logical Absolute
#include "real.fh"
*
      dR_A  = R_A-R
      dR_B  = R_B-R
      Call Multipole_E(q_A,Dipole_A,dR_A,E_A)
      Call Multipole_E(q_B,Dipole_B,dR_B,E_B)
*
      If (Absolute) Then
         Multipole_Expansion = Abs(E_A + E_B)
      Else
         Multipole_Expansion = E_A + E_B
      End If
      Return
      End
