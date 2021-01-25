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
      subroutine update_cvb(dx)
      implicit real*8 (a-h,o-z)
c ... Make: up to date? ...
      logical, external :: up2date_cvb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension dx(*)
c      dimension orbs(norb,norb),cvb(nvb)

      if(orbopt)call touch_cvb('ORBS')
      if(strucopt)call touch_cvb('CVB')
      call make_cvb('WFN')
c  TRY quantities only up2date if SVB/EVB ok :
      if(up2date_cvb('SVBTRY'))call make_cvb('SVB')
      if(up2date_cvb('EVBTRY'))call make_cvb('EVB')

      ic = 1
      i1 = mstackr_cvb(norb*norb)
      i2 = mstackr_cvb(nvb)
      i3 = mstackr_cvb(norb*norb)
      call update2_cvb(w(i1),w(i2),w(lv(1)),w(lv(2)),
     >  w(lw(2)),dx,
     >  ic,
     >  norb,nvb,nprorb,npr,orbopt,strucopt,sym,
     >  w(lp(6)),iw(ls(11)),nort,w(i3))
      call fmove_cvb(w(i1),w(lv(1)),norb*norb)
      call fmove_cvb(w(i2),w(lv(2)),nvb)
      call str2vbc_cvb(w(lv(2)),w(lv(5)))
      call mfreer_cvb(i1)
      return
      end
      subroutine upd_cvb(dx,orbs,cvb)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension dx(*)
      dimension orbs(norb,norb),cvb(nvb)
      if(orbopt)call touch_cvb('ORBSTRY')
      if(strucopt)call touch_cvb('CVBTRY')
      call make_cvb('WFNTRY')
      ic = 2
      i1 = mstackr_cvb(norb*norb)
      call update2_cvb(orbs,cvb,w(lv(1)),w(lv(2)),
     >  w(lw(2)),dx,
     >  ic,
     >  norb,nvb,nprorb,npr,orbopt,strucopt,sym,
     >  w(lp(6)),iw(ls(11)),nort,w(i1))
      call mfreer_cvb(i1)
      return
      end
