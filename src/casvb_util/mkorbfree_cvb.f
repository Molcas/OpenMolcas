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

#include "malloc_cvb.fh"

      i1 = mstackr_cvb(norb*norb)
      i2 = mstackr_cvb(norb*norb)
      i3 = mstackr_cvb(norb*norb)
      i4 = mstacki_cvb(nprorb)
      call mkorbfree2_cvb(w(lv(1)),iw(ls(3)),w(ls(4)),
     >  iw(ls(5)),w(ls(6)),iw(ls(8)),iw(ls(11)),iw(ls(12)),
     >  w(ls(14)),w(i1),w(i2),w(i3),iw(i4))
      call mfreer_cvb(i1)
      return
      end
