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
      Subroutine Get_D1MO(D1MO,nDens)
      Implicit Real*8 (A-H,O-Z)
#include "WrkSpc.fh"

      Character*24 Label
      Logical      Found
      Integer nDens
      Real*8 D1MO(nDens)

      Label='D1mo'
      Call qpg_dArray(Label,Found,mDens)
      If(.not.Found .or. nDens.eq.0) Then
         Call SysAbendMsg('get_d1mo','Did not find:',Label)
      End If
      If (mDens/=nDens) Then
         Write (6,*) 'Get_D1MO: mDens/=nDens'
         Write (6,*) 'mDens=',mDens
         Write (6,*) 'nDens=',nDens
         Call Abend()
      End If
      Call Get_dArray(Label,D1MO,nDens)

      Return
      End
