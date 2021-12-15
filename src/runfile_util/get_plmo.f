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
      Subroutine Get_PLMO(PLMO,nPLMO)
      Implicit Real*8 (A-H,O-Z)
      Character*24 Label
      Logical      Found
      Real*8 PLMO(nPLMO)

      Label='PLMO'
      Call qpg_dArray(Label,Found,mPLMO)
      If(.not.Found .or. mPLMO.eq.0) Then
         Call SysAbendMsg('get_plmo','Did not find:',Label)
      End If
      If (nPLMO/=mPLMO) Then
         Write (6,*) 'Get_PLMO: nPLMO/=mPLMO'
         Write (6,*) 'nPLMO=',nPLMO
         Write (6,*) 'mPLMO=',mPLMO
         Call Abend()
      End If
      Call get_dArray(Label,PLMO,nPLMO)

      Return
      End
