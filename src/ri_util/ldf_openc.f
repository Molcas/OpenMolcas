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
      Integer Function LDF_OpenC()
      Implicit None
      Integer LuC
#if defined (_MOLCAS_MPP_)
#include "para_info.fh"
      LuC=7
      If (nProcs.gt.1 .and. Is_Real_Par()) Then
         Call DAName_MF_WA(LuC,'MyLC')
      Else
         Call DAName_MF_WA(LuC,'LDFC')
      End If
#else
      LuC=7
      Call DAName_MF_WA(LuC,'LDFC')
#endif
      LDF_OpenC=LuC
      End
