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
      Subroutine Get_LCMO(LCMO,nLCMO)
      Implicit Real*8 (A-H,O-Z)
      Character*24 Label
      Logical      Found
      Real*8 LCMO(nLCMO)

      Label='LCMO'
      Call qpg_dArray(Label,Found,mLCMO)
      If(.not.Found .or. mLCMO.eq.0) Then
         Call SysAbendMsg('get_lcmo','Did not find:',Label)
      End If
      If (nLCMO/=mLCMO) Then
         Write (6,*) 'Get_LCMO: nLCMO/=mLCMO'
         Write (6,*) 'nLCMO=',nLCMO
         Write (6,*) 'mLCMO=',mLCMO
         Call Abend()
      End If
      Call get_dArray(Label,LCMO,nLCMO)

      Return
      End
