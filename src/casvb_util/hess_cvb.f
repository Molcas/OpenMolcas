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

#include "WrkSpc.fh"
      dimension vec(*)

      n_hess=n_hess+1
      if(.not.up2date_cvb('OOHESS'))then
        call make_cvb('OOHESS')
        call oohess_cvb(work(lv(1)),work(lc(2)),work(lc(3)),work(lc(4)),
     >    work(lw(1)),work(lw(2)),work(lw(3)),
     >    work(lq(8)),work(lq(9)),work(lq(3)),work(lq(4)))
      endif
      i1 = mstackr_cvb(npr)
      i2 = mstackr_cvb(npr)
      i3 = mstackr_cvb(norb*norb)
      i4 = mstackr_cvb(norb*norb)
      call free2all_cvb(vec,work(i1),1)
      if(icrit.eq.1)then
        call hess_svb1_cvb(work(lv(1)),
     >    work(lc(2)),work(lc(3)),work(lc(4)),work(lc(5)),
     >    work(lw(1)),work(lw(2)),work(lw(3)),
     >    work(lw(4)),work(lw(5)),work(lw(6)),
     >    work(lw(10)),
     >    work(lq(7)),work(lq(8)),work(lq(3)),
     >    work(lq(10)),iwork(ls(11)),
     >    work(i1),work(i2),work(i3),work(i4))
      elseif(icrit.eq.2)then
        call hess_evb1_cvb(work(lv(1)),
     >    work(lc(2)),work(lc(3)),work(lc(4)),
     >    work(lw(2)),work(lw(3)),
     >    work(lw(4)),work(lw(5)),work(lw(6)),
     >    work(lw(10)),
     >    work(lq(7)),work(lq(8)),work(lq(3)),
     >    work(lq(10)),iwork(ls(11)),
     >    work(i1),work(i2))
      endif
      call all2free_cvb(work(i2),vec,1)
      call mfreer_cvb(i1)
      return
      end
