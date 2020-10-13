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
      Subroutine Get_P2MO(P2MO,nP2MO)
      Implicit Real*8 (A-H,O-Z)
#include "WrkSpc.fh"

      Character*24 Label
      Logical      Found
      Integer nP2MO
      Real*8 P2MO(nP2MO)

      Label='P2MO'
      Call qpg_dArray(Label,Found,mP2MO)
      If(.not.Found .or. mP2MO.eq.0) Then
         Call SysAbendmsg('Get_P2MO','Did not find:',label)
      End If
      If (nP2MO/=mP2MO) Then
         Write (6,*) 'Get_P2MO: nP2MO/=mP2MO'
         Write (6,*) 'mP2MO=',mP2MO
         Write (6,*) 'nP2MO=',nP2MO
         Call Abend()
      End If
      Call get_dArray(Label,P2MO,nP2MO)

      Return
      End

      Subroutine Get_P2MOt(P2MO,nP2MO)
      Implicit Real*8 (A-H,O-Z)
#include "WrkSpc.fh"

      Character*24 Label
      Logical      Found
      Real*8 P2MO(nP2MO)

      Label='P2MOT'
      Call qpg_dArray(Label,Found,mP2MO)
      If(.not.Found .or. mP2MO.eq.0) Then
         Call SysAbendmsg('Get_P2MOt','Did not find:',label)
      End If
      If (nP2MO/=mP2MO) Then
         Write (6,*) 'Get_P2MO: nP2MO/=mP2MO'
         Write (6,*) 'mP2MO=',mP2MO
         Write (6,*) 'nP2MO=',nP2MO
         Call Abend()
      End If
      Call get_dArray(Label,P2MO,nP2MO)

      Return
      End
