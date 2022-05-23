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
      subroutine construc2_cvb(tconstr)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension tconstr(nvb,nvb)
      dimension dum(1)

      iconstruc_kp=iconstruc
      iconstruc=1
      irepm = mstackr_cvb(nvb)

      call span0_cvb(nvb,nvb)
      do 100 ivb=1,nvb
      call fzero(w(irepm),nvb)
      w(ivb+irepm-1)=-1d0
      call symtrizcvb_cvb(w(irepm))
      w(ivb+irepm-1)=w(ivb+irepm-1)+1d0
      call span1_cvb(w(irepm),1,dum,nvb,0)
100   continue
      call span2_cvb(tconstr,nconstr,dum,nvb,0)

      call mfreer_cvb(irepm)
      iconstruc=iconstruc_kp
      return
      end
