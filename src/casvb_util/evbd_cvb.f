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
      subroutine evbd_cvb(orbs,cvb,fx,ioptc,iter)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension orbs(norb,norb),cvb(nvb)

      mstackr_cvb1=mavailr_cvb()
      maxdav=min((mstackr_cvb1-ndet-nvb)/(3*nvb),nvb)
      i1 = mstackr_cvb(nvb*maxdav)
      i2 = mstackr_cvb(nvb*maxdav)
      i3 = mstackr_cvb(nvb*maxdav)
      i4 = mstackr_cvb(nvb)
      i5 = mstackr_cvb(maxdav*maxdav)
      i6 = mstackr_cvb(maxdav)
      i7 = mstackr_cvb(maxdav)
      call evbd2_cvb(orbs,cvb,fx,ioptc,iter,
     >  w(lw(4)),w(lw(5)),w(lw(6)),
     >  w(i1),w(i2),w(i3),w(i4),w(i5),w(i6),w(i7))
      call mfreer_cvb(i1)
      return
      end
