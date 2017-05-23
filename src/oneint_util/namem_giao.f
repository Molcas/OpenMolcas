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
      Subroutine NAMem_GIAO(nRys,MemNA_GIAO,la,lb,lr)
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
*
      labMin=nabSz(Max(la,lb)-1)+1
      labMax=nabSz(la+lb)
      lab = (labMax-labMin+1)
      kab = nElem(la)*nElem(lb)
*
      lcdMin_EF=nabSz(lr-1)+1
      lcdMax_EF=nabSz(lr)
      lcd_EF = (lcdMax_EF-lcdMin_EF+1)
      labcd_EF = lab*lcd_EF
      Mem0=labcd_EF
*
      lcdMin_NA=nabSz(lr-2)+1
      lcdMax_NA=nabSz(lr-1)
      lcd_NA = (lcdMax_NA-lcdMin_NA+1)
      labcd_NA = lab*lcd_NA
*
      Call mHRR(la,lb,nFlop,nMem)
      Mem1=Max(lcd_EF,lcd_NA)*nMem
*
      iAngV(1) = la
      iAngV(2) = lb
      iAngV(3) = lc
      iAngV(4) = ld
      Call MemRys(iAngV,Mem2_EF)
      iAngV(3) = 0
      Call MemRys(iAngV,Mem2_NA)
      Mem2 = Max(Mem2_EF,Mem2_NA,kab*Max(lcd_EF,lcd_NA))
*
      MemNA_GIAO = Mem0 + Mem1 + Mem2
      Return
      End
