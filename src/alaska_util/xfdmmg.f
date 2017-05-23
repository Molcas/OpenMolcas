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
      Subroutine XFdMmg(nRys,MnXFdG,la,lb,lr)
      Integer iAng(4)
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      MnXFdG=0
      Do iOrdOp = 0, 1
         iAng(1) = la
         iAng(2) = lb
         iAng(3) = iOrdOp
         iAng(4) = 0
         Call MemRg1(iAng,nRys,MemTmp)
         MemTmp = MemTmp + 2 + nElem(la)*nElem(lb)*nElem(iOrdOp)
         MnXFdG = Max(MnXFdG,MemTmp)
      End Do
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(lr)
      End
