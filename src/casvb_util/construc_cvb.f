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
c  *************************************************************
c  ** Routines for imposing constraints on VB wfn. parameters **
c  *************************************************************
c  *********************
c  ** Set-up routines **
c  *********************
      subroutine construc_cvb(tconstr,ipermzeta)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"

#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension tconstr(nvb,nvb),ipermzeta(norb,nzeta)

      i1 = mstackr_cvb(norb*norb)
      i2 = mstackr_cvb(norb*norb)
      i3 = mstackr_cvb(norb*norb)
      call setipermzeta_cvb(ipermzeta,
     >  w(lv(1)),w(ls(1)),iw(ls(13)),
     >  w(i1),w(i2),w(i3))
      call mfreer_cvb(i1)
      if(iconstruc.eq.2)call construc2_cvb(tconstr)
      return
      end
