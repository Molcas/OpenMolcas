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
      SubRoutine Cho_P_Distrib_Vec(Jin,Jfi,iDV,nV)
      Implicit none
      Integer  Jin, Jfi, nV
      Integer  iDV(*)
#include "cho_para_info.fh"

      Integer J0, J

      If (Cho_Real_Par) Then
         Call Cho_Distrib_Vec(Jin,Jfi,iDV,nV)
      Else
         J0 = Jin - 1
         nV = Jfi - J0
         Do J = 1,nV
            iDV(J) = J0 + J
         End Do
      End If

      End
