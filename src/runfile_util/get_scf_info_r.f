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
      Subroutine Get_SCF_Info_R(iCase,ipR,nR)
      Implicit Real*8 (A-H,O-Z)
#include "WrkSpc.fh"

      Character*24 Label
      Logical      Found

      Label='SCFInfoR'
      if(iCase.eq.1) Label='SCFInfoR_ab'
      Call qpg_dArray(label,Found,nR)
      If(.not.Found .or. nR.eq.0) then
         Call SysAbendmsg('get_scf_info_r','Did not find:',Label)
      End If
      Call GetMem('ipR','Allo','Real',ipR,nR)
      Call get_dArray(label,Work(ipR),nR)

      Return
      End
