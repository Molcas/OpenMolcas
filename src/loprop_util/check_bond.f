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
      Logical Function Check_Bond(CXi,CXj,iANr,jANr,Factor)
* Returns true if the bond should be included, otherwise false
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 CXi(3),CXj(3)
*
      Check_Bond = .True.
      If (Factor .lt. Zero) Then
         Check_Bond = .True.
      Else
         Radius_i    = Bragg_Slater(iANr)
         Radius_j    = Bragg_Slater(jANr)
         Bond_Length = Sqrt((CXi(1)-CXj(1))**2
     &                     +(CXi(2)-CXj(2))**2
     &                     +(CXi(3)-CXj(3))**2)
         Bond_Max    = Factor*(Radius_i+Radius_j)
         If (Bond_Length .gt. Bond_Max) Check_Bond = .False.
      End If
*
      Return
      End
