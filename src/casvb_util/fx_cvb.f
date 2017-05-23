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
      subroutine fx_cvb(fx,fast)
c
c  FAST=.TRUE. :     Quick evaluation of f(x) --
c                    no reusable quantities
c  FAST=.FALSE.:     Slower evaluation of f(x) --
c                    but with reusable quantities
c                    (CIVECP/CIVBH & CIVBS)
c
      implicit real*8 (a-h,o-z)
      logical fast
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"

      dxmove=.true.
      iwf1 = lv(1)
      iwf2 = lv(2)
      if(.not.memplenty)then
        call ciwr_cvb(w(lc(2)),61002.2d0)
        call ciwr_cvb(w(lc(3)),61003.2d0)
        call ciwr_cvb(w(lc(4)),61004.2d0)
        call setcnt2_cvb(2,0)
        call setcnt2_cvb(3,0)
        call setcnt2_cvb(4,0)
      endif
      call setcnt2_cvb(6,0)
      call setcnt2_cvb(7,0)
      call setcnt2_cvb(8,0)
      if(icrit.eq.1)then
        call fx_svb1_cvb(fx,fast,w(iwf1),w(iwf2),
     >    w(lc(1)),w(lc(6)),w(lc(7)),w(lc(8)),
     >    w(lw(4)),w(lw(5)),w(lw(6)),w(lw(9)))
      elseif(icrit.eq.2)then
        call fx_evb1_cvb(fx,fast,w(iwf1),w(iwf2),
     >    w(lc(1)),w(lc(6)),w(lc(7)),w(lc(8)),
     >    w(lw(4)),w(lw(5)),w(lw(6)),w(lw(9)))
      endif
      if(.not.memplenty)then
        call ciwr_cvb(w(lc(6)),61006.2d0)
        call ciwr_cvb(w(lc(7)),61007.2d0)
        call ciwr_cvb(w(lc(8)),61008.2d0)
        call cird_cvb(w(lc(2)),61002.2d0)
        call cird_cvb(w(lc(3)),61003.2d0)
        call cird_cvb(w(lc(4)),61004.2d0)
      endif
c  Figure out what we just calculated, and make it up2date :
      if(dxmove)then
        if(icrit.eq.1)then
          call make_cvb('SVB')
        elseif(icrit.eq.2)then
          call make_cvb('EVB')
        endif
      else
        if(icrit.eq.1)then
          call make_cvb('SVBTRY')
        elseif(icrit.eq.2)then
          call make_cvb('EVBTRY')
        endif
      endif
      return
      end
      subroutine fxdx_cvb(fx,fast,dx)
      implicit real*8 (a-h,o-z)
      logical fast
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension dx(*)

      dxmove=.false.
      iwf1 = lw(12)
      iwf2 = lw(13)
      call upd_cvb(dx,w(iwf1),w(iwf2))
      if(.not.memplenty)then
        call ciwr_cvb(w(lc(2)),61002.2d0)
        call ciwr_cvb(w(lc(3)),61003.2d0)
        call ciwr_cvb(w(lc(4)),61004.2d0)
        call setcnt2_cvb(2,0)
        call setcnt2_cvb(3,0)
        call setcnt2_cvb(4,0)
      endif
      call setcnt2_cvb(6,0)
      call setcnt2_cvb(7,0)
      call setcnt2_cvb(8,0)
      if(icrit.eq.1)then
        call fx_svb1_cvb(fx,fast,w(iwf1),w(iwf2),
     >    w(lc(1)),w(lc(6)),w(lc(7)),w(lc(8)),
     >    w(lw(4)),w(lw(5)),w(lw(6)),w(lw(9)))
      elseif(icrit.eq.2)then
        call fx_evb1_cvb(fx,fast,w(iwf1),w(iwf2),
     >    w(lc(1)),w(lc(6)),w(lc(7)),w(lc(8)),
     >    w(lw(4)),w(lw(5)),w(lw(6)),w(lw(9)))
      endif
      if(.not.memplenty)then
        call ciwr_cvb(w(lc(6)),61006.2d0)
        call ciwr_cvb(w(lc(7)),61007.2d0)
        call ciwr_cvb(w(lc(8)),61008.2d0)
        call cird_cvb(w(lc(2)),61002.2d0)
        call cird_cvb(w(lc(3)),61003.2d0)
        call cird_cvb(w(lc(4)),61004.2d0)
      endif
c  Figure out what we just calculated, and make it up2date :
      if(dxmove)then
        if(icrit.eq.1)then
          call make_cvb('SVB')
        elseif(icrit.eq.2)then
          call make_cvb('EVB')
        endif
      else
        if(icrit.eq.1)then
          call make_cvb('SVBTRY')
        elseif(icrit.eq.2)then
          call make_cvb('EVBTRY')
        endif
      endif
      return
      end
