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
      subroutine mksyminit_cvb()
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"

c  Generate symmetry information - first call gets dimensions :
      i1 = mstacki_cvb(4*norbrel)
      i2 = mstacki_cvb(norb*norbrel)
      i3 = mstacki_cvb(norb)
      i4 = mstackr_cvb(norb*norb)
      i5 = mstackr_cvb(norb*norb)
      i6 = mstackr_cvb(norb)
      i7 = mstackr_cvb(norb)
      i8 = mstackr_cvb(norb*norb)
      i9 = mstackr_cvb(norb*norb)
      i10= mstacki_cvb(norb)
      call syminit2_cvb(w(ls(1)),iw(ls(2)),iw(ls(3)),w(ls(4)),
     >  iw(ls(5)),w(ls(6)),iw(i1),iw(i2),iw(i3),w(i4),w(i5),
     >  w(i6),w(i7),w(i8),w(i9),iw(i10),
     >  iw(ls(8)),iw(ls(9)),iw(ls(10)),
     >  iw(ls(11)),iw(ls(12)),iw(ls(13)))
      call mfreei_cvb(i1)
      return
      end
