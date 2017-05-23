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
      subroutine symtrizorbs_cvb(orbs)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension orbs(norb,norb)

      i1 = mstacki_cvb(norb)
      i2 = mstacki_cvb(norb)
      i3 = mstacki_cvb(norb)
      i4 = mstacki_cvb(norb)
      i5 = mstacki_cvb(norb)
      i6 = mstackr_cvb(norb)
      i7 = mstackr_cvb(norb)
      call symtrizorbs2_cvb(orbs,
     >  iw(ls(3)),w(ls(4)),iw(ls(5)),w(ls(6)),iw(ls(8)),iw(ls(11)),
     >  iw(i1),iw(i2),iw(i3),iw(i4),iw(i5),w(i6),w(i7))
      call mfreei_cvb(i1)
      return
      end
