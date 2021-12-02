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
      subroutine bspset_cvb(kbasis1,ic,need)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"


#include "WrkSpc.fh"

      if(ic.eq.1)then
        i1 = mstackiz_cvb((nel+1)*(nel+1)*(nel+1))
        call bspset2_cvb(iwork(i1),nel,kbasis1,need)
        call mfreei_cvb(i1)
      elseif(ic.eq.2)then
        do 100 i=0,(nel+1)*(nel+1)*(nel+1)-1
        iwork(i+lb(3))=-1
100     continue
        call bspset2_cvb(iwork(lb(3)),nel,kbasis1,need)
        call setifnss_cvb(iwork(lb(4)),iwork(lb(5)),iwork(lb(6)))
      endif
      if(kbasis1.eq.6)need=0
      return
      end
