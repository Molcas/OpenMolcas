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
* Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
*               1996-2006, David L. Cooper                             *
************************************************************************
      subroutine sorti_cvb(n,arrin)
      implicit real*8 (a-h,o-z)
#include "malloc_cvb.fh"
      integer arrin
      dimension arrin(n)

      i1 = mstacki_cvb(n)
      call sortindxi_cvb(n,arrin,iw(i1))
      i2 = mstacki_cvb(n)
      do 100 i=0,n-1
      iw(i+i2)=arrin(iw(i+i1))
100   continue
      call imove_cvb(iw(i2),arrin,n)
      call mfreei_cvb(i1)
      return
      end
