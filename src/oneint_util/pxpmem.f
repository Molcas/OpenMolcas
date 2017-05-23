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
      Subroutine PXPMem(nRys,MemPXP,la,lb,lr)
*
*     Statement function for Cartesian index
*
      MemPXP=0
      nRys =0
      Call PXMem(nOrder,Mem,la,lb+1,lr-1)
      MemPXP=Max(MemPXP,Mem)
      nRys =Max(nRys,nOrder)
      If (lb.GT.0) Then
         Call PXMem(nOrder,Mem,la,lb-1,lr-1)
         MemPXP=Max(MemPXP,Mem)
         nRys =Max(nRys,nOrder)
      End If
*
      Return
      End
