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
* Copyright (C) 1991, Roland Lindh                                     *
************************************************************************
      Subroutine PVMem(nRys,MemPV,la,lb,lr,KrnMem)
      External KrnMem
*
      Call KrnMem(nRys,MemNA1,la+1,lb,lr-1)
*
      If (la.ne.0) Then
         Call KrnMem(nRys,MemNA2,la-1,lb,lr-1)
      Else
         MemNA2=0
      End If
*
      MemPV=Max(MemNA1,MemNA2)
*
      Return
      End
