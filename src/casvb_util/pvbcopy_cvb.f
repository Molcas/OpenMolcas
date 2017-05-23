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
      subroutine pvbcopy_cvb(cfrom,cto)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension cfrom(*),cto(*)

      icfrom=nint(cfrom(1))
      icto=nint(cto(1))
      if(iform_ci(icfrom).ne.0.or.iform_ci(icto).ne.0)then
        write(6,*)' Unsupported format in PVBCOPY'
        call abend_cvb()
      endif
      call pvbcopy2_cvb(w(iaddr_ci(icfrom)),w(iaddr_ci(icto)),
     >  iw(ll(11)),iw(ll(12)),dum,0)
      call setcnt2_cvb(icto,0)
      return
      end
      subroutine pvbdot_cvb(cfrom,cto,ret)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension cfrom(*),cto(*)
      icfrom=nint(cfrom(1))
      icto=nint(cto(1))
      if(iform_ci(icfrom).ne.0.or.iform_ci(icto).ne.0)then
        write(6,*)' Unsupported format in PVBDOT'
        call abend_cvb()
      endif
      call pvbcopy2_cvb(w(iaddr_ci(icfrom)),w(iaddr_ci(icto)),
     >  iw(ll(11)),iw(ll(12)),ret,1)
      return
      end
