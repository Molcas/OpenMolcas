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
      Subroutine Get_DLMO(ipDLMO,nDens)
      Implicit Real*8 (A-H,O-Z)
#include "WrkSpc.fh"
      Character*24 Label
      Logical      Found

      Label='DLMO'
      Call qpg_dArray(Label,Found,nDens)
      If(.not.Found .or. nDens.eq.0) Then
         Call SysAbendMsg('get_dlmo','Did not find:',Label)
      End If
      Call GetMem('Dens','Allo','Real',ipDLMO,nDens)
      Call Get_dArray(Label,Work(ipDLMO),nDens)

      Return
      End
