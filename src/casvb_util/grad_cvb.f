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
      subroutine grad_cvb(grad)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "fx_cvb.fh"
#include "malloc_cvb.fh"

      dimension grad(*)

      call touch_cvb('OOHESS')
      if(dxmove.and.memplenty)then
        call cicopy_cvb(w(lc(6)),w(lc(2)))
        call cicopy_cvb(w(lc(7)),w(lc(3)))
        call cicopy_cvb(w(lc(8)),w(lc(4)))
      elseif(dxmove)then
        call cird_cvb(w(lc(2)),61006.2d0)
        call cird_cvb(w(lc(3)),61007.2d0)
        call cird_cvb(w(lc(4)),61008.2d0)
      endif
      ovraa=ovraa_try
      ovrab=ovrab_try
      ww=ww_try
      if(icrit.eq.1)then
        call gr_svb1_cvb(w(lc(2)),w(lc(3)),w(lc(4)),w(lw(10)),
     >    grad,w(lq(7)),w(lq(8)),w(lq(9)),w(lq(10)))
      elseif(icrit.eq.2)then
        call gr_evb1_cvb(w(lc(2)),w(lc(3)),w(lc(4)),w(lw(10)),
     >    grad,w(lq(7)),w(lq(8)),w(lq(9)),w(lq(10)))
      endif
      return
      end
