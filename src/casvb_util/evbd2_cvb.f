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
      subroutine evbd2_cvb(orbs,cvb,fx,ioptc,iter,
     >  gjorb,gjorb2,gjorb3,
     >  c,axc,sxc,res,hp,solp,solp_res)
      implicit real*8 (a-h,o-z)
      external asonc_cvb,ddsol7_cvb,ddres7_cvb,ddres2upd10_cvb
      external ddrestart_cvb
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "davtune_cvb.fh"
#include "opt2_cvb.fh"
#include "malloc_cvb.fh"
      dimension orbs(norb,norb),cvb(nvb)
      dimension gjorb(*),gjorb2(*),gjorb3(*)
      dimension c(nvb,maxdav),axc(nvb,maxdav),sxc(nvb,maxdav),res(nvb)
      dimension hp(maxdav,maxdav),solp(maxdav),solp_res(maxdav)
      dimension dum(max(nvb,maxdav))

      call makegjorbs_cvb(orbs,gjorb,gjorb2,gjorb3)

      ioptc=1
      nvguess=1
      nvrestart=0
      call fmove_cvb(cvb,c,nvb)
      if(.not.follow)then
        ifollow=2
      else
        ifollow=4
      endif

      call ddinit7_cvb(ifollow,isaddle,ip(3))
      call ddres2updinit_cvb(0)
      call dirdiag_cvb(
     >  asonc_cvb,ddsol7_cvb,ddres7_cvb,ddres2upd10_cvb,
     >  ddrestart_cvb,
     >  c,axc,sxc,.false.,cvb,res,dum,
     >  hp,dum,solp,solp_res,
     >  .false.,.true.,.false.,maxdav,nvb,nvb,
     >  nvguess,nvrestart,isaddle,ifollow,mxiter,
     >  resthr,orththr,nortiter,corenrg,
     >  ioptc,iter,fx,ip(3))
      have_solved_it=.true.

      ovraa=one
      evb=fx
      return
      end
