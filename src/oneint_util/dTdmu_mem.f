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
      Subroutine dTdmu_mem(nRys,Mem_dTdMu,la,lb,lr)
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      Mem_dTdMu=0
      nRys =0
      Call EFMmP(nOrder,MmEFP,la,lb+1,lr)
      Mem_dTdMu=Max(Mem_dTdMu,MmEFP)
      nRys =Max(nRys,nOrder)
      If (lb.ge.1) Then
         Call EFMmP(nOrder,MmEFP,la,lb-1,lr)
         Mem_dTdMu=Max(Mem_dTdMu,MmEFP)
         nRys =Max(nRys,nOrder)
      End If
*
*     Add a scratch area for intermediate integrals
*
      MemDer = 3*nElem(la)*nElem(lb+1)
      If (lb.ge.1) MemDer=MemDer + 3*nElem(la)*nElem(lb-1)
      Mem_dTdMu = Mem_dTdMu + MemDer + 1
      Mem_dTdMu = Mem_dTdMu + nElem(la)*nElem(lb)*3
*
      Return
      End
