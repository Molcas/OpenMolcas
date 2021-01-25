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
      subroutine hess_cvb(vec)
      implicit real*8 (a-h,o-z)
c ... Make: up to date? ...
      logical, external :: up2date_cvb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension vec(*)

      n_hess=n_hess+1
      if(.not.up2date_cvb('OOHESS'))then
        call make_cvb('OOHESS')
        call oohess_cvb(w(lv(1)),w(lc(2)),w(lc(3)),w(lc(4)),
     >    w(lw(1)),w(lw(2)),w(lw(3)),
     >    w(lq(8)),w(lq(9)),w(lq(3)),w(lq(4)))
      endif
      i1 = mstackr_cvb(npr)
      i2 = mstackr_cvb(npr)
      i3 = mstackr_cvb(norb*norb)
      i4 = mstackr_cvb(norb*norb)
      call free2all_cvb(vec,w(i1),1)
      if(icrit.eq.1)then
        call hess_svb1_cvb(w(lv(1)),
     >    w(lc(2)),w(lc(3)),w(lc(4)),w(lc(5)),
     >    w(lw(1)),w(lw(2)),w(lw(3)),
     >    w(lw(4)),w(lw(5)),w(lw(6)),
     >    w(lw(10)),
     >    w(lq(7)),w(lq(8)),w(lq(3)),
     >    w(lq(10)),iw(ls(11)),
     >    w(i1),w(i2),w(i3),w(i4))
      elseif(icrit.eq.2)then
        call hess_evb1_cvb(w(lv(1)),
     >    w(lc(2)),w(lc(3)),w(lc(4)),
     >    w(lw(2)),w(lw(3)),
     >    w(lw(4)),w(lw(5)),w(lw(6)),
     >    w(lw(10)),
     >    w(lq(7)),w(lq(8)),w(lq(3)),
     >    w(lq(10)),iw(ls(11)),
     >    w(i1),w(i2))
      endif
      call all2free_cvb(w(i2),vec,1)
      call mfreer_cvb(i1)
      return
      end
