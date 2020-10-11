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
      Subroutine Get_D1AV(D1AV,nD1AV)
      Implicit Real*8 (A-H,O-Z)
      Character*24 Label
      Logical      Found
      Real*8 D1AV(nD1AV)

      Label='D1av'
      Call qpg_dArray(Label,Found,mD1AV)
      If(.not.Found .or. mD1AV.eq.0) Then
         Call SysAbendMsg('Get_D1AV','Did not find:',Label)
      End If
      If (nD1AV/=mD1AV) Then
         Write (6,*) 'Get_D1AV: nD1AV/=mD1AV'
         Write (6,*) 'nD1AV=',nD1AV
         Write (6,*) 'mD1AV=',mD1AV
         Call Abend()
      End If

      Call Get_dArray(Label,D1AV,nD1AV)

      Return
      End
