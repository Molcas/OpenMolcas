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
      Subroutine EFMmP(nRys,MmEFP,la,lb,lr)
*
      Integer iAngV(4)
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6  - 1
*
      lc = lr
      ld = 0
      nRys = (la+lb+lc+ld+2)/2
      labMin=nabSz(Max(la,lb)-1)+1
      labMax=nabSz(la+lb)
      lcdMin=nabSz(lr-1)+1
      lcdMax=nabSz(lr)
      lab = (labMax-labMin+1)
      kab = nElem(la)*nElem(lb)
      lcd = (lcdMax-lcdMin+1)
      labcd = lab*lcd
*
      Call mHRR(la,lb,nFlop,nMem)
      Mem1=Max(lcd*nMem,labcd)
*
      iAngV(1) = la
      iAngV(2) = lb
      iAngV(3) = lc
      iAngV(4) = ld
      Call MemRys(iAngV,Mem2)
      Mem2 = Max(Mem2,kab*lcd)
*
      MmEFP=Mem1 + Mem2
      Return
      End
