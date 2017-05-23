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
      Subroutine Phi_point(iPhi,nPhi,Cos_Phi,Sin_Phi,w_Phi)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
*
      q = Pi*(Two*DBLE(iPhi)-1.0d0)/DBLE(nPhi)
      If (Abs(Cos(q)).gt.1.0D-14) Then
         Cos_Phi=Cos(q)
      Else
         Cos_Phi=Zero
      End If
      If (Abs(Sin(q)).gt.1.0D-14) Then
         Sin_Phi=Sin(q)
      Else
         Sin_Phi=Zero
      End If
      w_Phi=Two*Pi/DBLE(nPhi)
*
      Return
      End
