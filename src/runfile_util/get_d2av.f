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
      Subroutine Get_D2AV(D2AV,nD2AV)
      Implicit Real*8 (A-H,O-Z)
      Character*24 Label
      Logical      Found
      Real*8 D2AV(nD2AV)

      Label='D2av'
      Call qpg_dArray(Label,Found,mD2AV)
      If(.not.Found .or. mD2AV.eq.0) Then
         Call SysAbendMsg('get_d2av','Did not find',Label)
      End If
      If (nD2AV/=mD2AV) Then
         Write (6,*) 'Get_D2AV: nD2AV/=mD2AV'
         Write (6,*) 'nD2AV=',nD2AV
         Write (6,*) 'mD2AV=',mD2AV
         Call Abend()
      End If
      Call get_dArray(Label,D2AV,nD2AV)

      Return
      End
