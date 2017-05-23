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
      subroutine asonc12e_cvb(c,axc,sxc,nvec,nprm)
c  Applies S and H on c vector(s).
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension c(nprm,nvec),axc(nprm,nvec),sxc(nprm,nvec)

      i1 = mstackr_cvb(nvb+nprorb)
      call asonc12e2_cvb(c,axc,sxc,nvec,nprm,
     >  w(lc(3)),w(lc(4)),w(lc(2)),
     >  w(lv(1)),w(lw(4)),w(lw(5)),w(lw(6)),w(lw(9)),
     >  w(lv(2)),
     >  w(i1))
      call mfreer_cvb(i1)
      return
      end
