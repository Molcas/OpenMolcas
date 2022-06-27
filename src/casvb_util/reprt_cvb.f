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
      subroutine reprt_cvb()
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "WrkSpc.fh"

      i1 = mstackr_cvb(norb*norb)
      i2 = mstackr_cvb(norb)

      if(lciweights) then
        lc_citmp=lc(5)
      else
C citmp is not used
        lc_citmp=lc(1)
      endif

      call reprt2_cvb(work(lv(1)),work(lv(2)),
     >  work(lc(1)),work(lc(2)),work(lc(3)),work(lc(4)),work(lc_citmp),
     >  work(lp(1)),work(lp(2)),
     >  work(lw(1)),work(lw(2)),work(lw(3)),work(lw(4)),work(lw(5)),
     >  work(lw(6)),
     >  work(lw(7)),work(lw(8)),work(lw(9)),work(lw(10)),work(lw(11)),
     >  work(i1),work(i2))
      call mfreer_cvb(i1)
      return
      end
