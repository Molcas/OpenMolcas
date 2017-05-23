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
* Copyright (C) 1996-2006, T. Thorsteinsson and D. L. Cooper           *
************************************************************************
      subroutine stringen_cvb(norb,nel,locc,lunocc,nstring)
      implicit real*8 (a-h,o-z)
#include "malloc_cvb.fh"
      dimension locc(*),lunocc(*)

      i_nkmin = mstacki_cvb(norb+1)
      i_nkmax = mstacki_cvb(norb+1)
c Spin string loop initialization (use xdet as graph storage) :
      do 100 iorb=0,norb
      iw(iorb+i_nkmin)=max(iorb-norb+nel,0)
100   iw(iorb+i_nkmax)=min(iorb,nel)
      call mmstringen_cvb(norb,nel,locc,lunocc,nstring,
     >  iw(i_nkmin),iw(i_nkmax))
      call mfreei_cvb(i_nkmin)
      return
      end
