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
      Subroutine Get_DLAO(DLAO,nDLAO)
      Implicit Real*8 (A-H,O-Z)
      Character*24 Label
      Logical      Found
      Real*8 DLAO(nDLAO)

      Label='DLAO'
      Call qpg_dArray(Label,Found,mDLAO)
      If(.not.Found .or. mDLAO.eq.0) Then
         Call SysAbendMsg('get_dlao','Did not find:',Label)
      End If
      If (nDLAO/=mDLAO) Then
         Write (6,*) 'Get_DLAO: nDLAO/=mDLAO'
         Write (6,*) 'nDLAO=',DLAOV
         Write (6,*) 'mDLAO=',mDLAO
         Call Abend()
      End If

      Call Get_dArray(Label,DLAO,nDLAO)

      Return
      End
