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
      Subroutine QpVMem(nHer,MemQpV,la,lb,lr)
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      Call MltMmP(nOrder,MemQpV,la,lb+1,lr-1)
      nHer = nOrder
      If (lb.gt.0) Then
         Call MltMmP(nOrder,MmMltP,la,lb-1,lr-1)
         MemQpV = Max(MemQpV,MmMltP) + nElem(la)*nElem(lb-1)*3
      End If
      MemQpV = MemQpV + 1 + nElem(la)*nElem(lb+1)*3
      MemQpV = MemQpV + nElem(la)*nElem(lb)*6
*
      Return
      End
