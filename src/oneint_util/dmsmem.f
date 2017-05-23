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
      Subroutine DMSMem(nRys,MemDMS,la,lb,lr)
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      MemDMS=0
      nRys =0
      Call EFMmP(nOrder,MmEFP,la,lb+1,lr-1)
      MemDMS=Max(MemDMS,MmEFP)
      nRys =Max(nRys,nOrder)
      Call EFMmP(nOrder,MmEFP,la,lb,lr-1)
      MemDMS=Max(MemDMS,MmEFP)
      nRys =Max(nRys,nOrder)
*
*     Add a scratch area for intermediate integrals
*
      MemDer = 3*(nElem(la)*nElem(lb+1) +  nElem(la)*nElem(lb))
      MemDMS = MemDMS + MemDer
      MemDMS = MemDMS + nElem(la)*nElem(lb)*9
*
      Return
      End
