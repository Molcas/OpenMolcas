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
      Subroutine DMSMem(
#define _CALLING_
#include "mem_interface.fh"
     &)
#include "mem_interface.fh"
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      Mem=0
      nHer =0
      Call EFMmP(nOrder,MmEFP,la,lb+1,lr-1)
      Mem=Max(Mem,MmEFP)
      nHer =Max(nHer,nOrder)
      Call EFMmP(nOrder,MmEFP,la,lb,lr-1)
      Mem=Max(Mem,MmEFP)
      nHer =Max(nHer,nOrder)
*
*     Add a scratch area for intermediate integrals
*
      MemDer = 3*(nElem(la)*nElem(lb+1) +  nElem(la)*nElem(lb))
      Mem = Mem + MemDer
      Mem = Mem + nElem(la)*nElem(lb)*9
*
      Return
      End
