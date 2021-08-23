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
      subroutine grad_cvb(grad)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "fx_cvb.fh"
#include "WrkSpc.fh"

      dimension grad(*)

      call touch_cvb('OOHESS')
      if(dxmove.and.memplenty)then
        call cicopy_cvb(work(lc(6)),work(lc(2)))
        call cicopy_cvb(work(lc(7)),work(lc(3)))
        call cicopy_cvb(work(lc(8)),work(lc(4)))
      elseif(dxmove)then
        call cird_cvb(work(lc(2)),61006.2d0)
        call cird_cvb(work(lc(3)),61007.2d0)
        call cird_cvb(work(lc(4)),61008.2d0)
      endif
      ovraa=ovraa_try
      ovrab=ovrab_try
      ww=ww_try
      if(icrit.eq.1)then
        call gr_svb1_cvb(work(lc(2)),work(lc(3)),work(lc(4)),
     >    work(lw(10)),
     >    grad,work(lq(7)),work(lq(8)),work(lq(9)),work(lq(10)))
      elseif(icrit.eq.2)then
        call gr_evb1_cvb(work(lc(2)),work(lc(3)),work(lc(4)),
     >    work(lw(10)),
     >    grad,work(lq(7)),work(lq(8)),work(lq(9)),work(lq(10)))
      endif
      return
      end
