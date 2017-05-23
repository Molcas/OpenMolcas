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
      Subroutine NAMmG(nRys,MmNAG,la,lb,lr)
      Integer iAng(4)
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      iAng(1) = la
      iAng(2) = lb
      iAng(3) = 0
      iAng(4) = 0
      Call MemRg1(iAng,nRys,MmNAG)
      MmNAG = MmNAG + 2 + nElem(la)*nElem(lb)
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(lr)
      End
