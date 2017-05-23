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
      subroutine reprt_cvb()
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"

      i1 = mstackr_cvb(norb*norb)
      i2 = mstackr_cvb(norb)

      if(lciweights) then
        lc_citmp=lc(5)
      else
C citmp is not used
        lc_citmp=lc(1)
      endif

      call reprt2_cvb(w(lv(1)),w(lv(2)),
     >  w(lc(1)),w(lc(2)),w(lc(3)),w(lc(4)),w(lc_citmp),
     >  w(lp(1)),w(lp(2)),
     >  w(lw(1)),w(lw(2)),w(lw(3)),w(lw(4)),w(lw(5)),w(lw(6)),
     >  w(lw(7)),w(lw(8)),w(lw(9)),w(lw(10)),w(lw(11)),
     >  w(i1),w(i2))
      call mfreer_cvb(i1)
      return
      end
