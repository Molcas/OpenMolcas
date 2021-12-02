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
      subroutine mkorbfree_cvb()
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "WrkSpc.fh"

      i1 = mstackr_cvb(norb*norb)
      i2 = mstackr_cvb(norb*norb)
      i3 = mstackr_cvb(norb*norb)
      i4 = mstacki_cvb(nprorb)
      call mkorbfree2_cvb(work(lv(1)),iwork(ls(3)),work(ls(4)),
     >  iwork(ls(5)),work(ls(6)),iwork(ls(8)),iwork(ls(11)),
     >  iwork(ls(12)),
     >  work(ls(14)),work(i1),work(i2),work(i3),iwork(i4))
      call mfreer_cvb(i1)
      return
      end
