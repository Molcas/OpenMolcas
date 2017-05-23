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
      Subroutine Get_SCF_Info_I(iCase,ipI,nI)
      Implicit Real*8 (A-H,O-Z)
#include "WrkSpc.fh"

      Character*24 Label
      Logical      Found

      Label='SCFInfoI'
      if(iCase.eq.1) Label='SCFInfoI_ab'
      Call qpg_iArray(Label,Found,nI)
      If(.not.Found .or. nI.eq.0) Then
         Call SysAbendMsg('get_scf_info_i','Did not find:',Label)
      End If
      Call GetMem('ipI','Allo','Inte',ipI,nI)
      Call Get_iArray(Label,iWork(ipI),nI)

      Return
      End
