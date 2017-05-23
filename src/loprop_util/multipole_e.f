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
      Subroutine Multipole_E(q,Dipole,R,E)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
*
      R_Inv    = 1.0D0/R
      R2_Inv   = R_Inv*R_Inv
*
      If (R .lt. Zero) Then
         xMonopole = -q*R2_Inv
         xDipole   = -3.0D0*Dipole*R*R2_Inv*R2_Inv
      Else
         xMonopole = q*R2_Inv
         xDipole   = 3.0D0*Dipole*R*R2_Inv*R2_Inv
      End If
*
      E = xMonopole + xDipole
*
      Return
      End
