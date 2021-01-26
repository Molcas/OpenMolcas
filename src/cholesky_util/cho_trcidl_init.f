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
      Use Para_Info, Only: nProcs
      use ChoArr, only: Idle
      Implicit None
#include "cho_para_info.fh"
#include "stdalloc.fh"

      Integer l_Idle

      If (Cho_Real_Par) Then
         l_Idle=nProcs
      Else
         l_Idle=1
      End If
      Call mma_allocate(Idle,l_Idle,Label='Idle')

      Idle(:)=0

      End
