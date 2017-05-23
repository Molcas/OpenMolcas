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
      SubRoutine Cho_TrcIdl_Init()
C
C     Thomas Bondo Pedersen, May 2010.
C
C     Allocate and init array for tracing idle processors
C
      Implicit None
#include "choptr2.fh"
#include "cho_para_info.fh"
#include "WrkSpc.fh"

      Integer i

      If (Cho_Real_Par) Then
         l_Idle=nProcs
      Else
         l_Idle=1
      End If
      Call GetMem('TrcIdl','Allo','Inte',ip_Idle,l_Idle)

      Do i=0,l_Idle-1
         iWork(ip_Idle+i)=0
      End Do

      End
