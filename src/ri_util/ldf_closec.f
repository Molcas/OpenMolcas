************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      Integer Function LDF_CloseC(LuC)
      Implicit None
      Integer LuC
#include "para_info.fh"
      If (LuC.lt.1) Then
         LDF_CloseC=-1
         Return
      End If
#if defined (_MOLCAS_MPP_)
      If (nProcs.gt.1 .and. Is_Real_Par()) Then
         Call DAEras(LuC)
      Else
         Call DAClos(LuC)
      End If
#else
      Call DAClos(LuC)
#endif
      LDF_CloseC=0
      End
