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
* Copyright (C) 1990, IBM                                              *
************************************************************************
      Subroutine  mHrr(la,lb,nSize,nMem)
      Implicit Real*8 (a-h,o-z)
*
#include "print.fh"
#include "WrkSpc.fh"
*
*     Statement function for canonical indices
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      iRout =25
      iPrint = nPrint(iRout)
*
*     First find the size of the working array.
*
      nMem  = 0
      nMem2 = 0
      nSize = 0
      Do 100 ib = 0, Min(la,lb)
         nMem1 = 0
         Do 150 ia = Max(la,lb), la+lb-ib
            nSize = nSize + nElem(ib)*nElem(ia)
            nMem1 = nMem1 + nElem(ib)*nElem(ia)
 150     Continue
         nMem=Max(nMem,nMem1+nMem2)
         nMem2 = nMem1
         If (ib.eq.0) nSize = 0
 100  Continue
*
      Return
      End
