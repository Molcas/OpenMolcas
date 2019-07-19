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
      subroutine svbd2_cvb(orbs,cvb,fx,ioptc,iter,
     >  civec,civbs,gjorb,gjorb2,gjorb3,cvbdet,
     >  c,sxc,res,rhs,rhsp,solp,solp_res)
      implicit real*8 (a-h,o-z)
      external asonc1_cvb,ddsolsvb_cvb,ddressvb_cvb,
     >  ddres2upd10_cvb
      external ddrestart_cvb
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "davtune_cvb.fh"
#include "opt2_cvb.fh"
      dimension orbs(norb,norb),cvb(nvb)
      dimension civec(ndet),civbs(ndet)
      dimension gjorb(*),gjorb2(*),gjorb3(*)
      dimension cvbdet(ndetvb)
      dimension c(nvb,maxdav),sxc(nvb,maxdav),res(nvb),rhs(nvb)
      dimension rhsp(maxdav),solp(maxdav),solp_res(maxdav)
      dimension dum(max(nvb,maxdav),maxdav)

      call makegjorbs_cvb(orbs,gjorb,gjorb2,gjorb3)

      if(memplenty)then
        call cicopy_cvb(civec,civbs)
      else
        call cird_cvb(civbs,61001.2d0)
      endif
      call applyt_cvb(civbs,gjorb2)
      call ci2vbg_cvb(civbs,cvbdet)
      call vb2strg_cvb(cvbdet,rhs)

      ioptc=1
      nvguess=1
      nvrestart=0
      call fmove_cvb(cvb,c,nvb)
      if(.not.follow)then
        ifollow=1
      else
        ifollow=4
      endif

      call ddinitsvb_cvb(ifollow,isaddle,ip(3))
      call ddres2updinit_cvb(0)
      call dirdiag_cvb(
     >  asonc1_cvb,ddsolsvb_cvb,ddressvb_cvb,ddres2upd10_cvb,
     >  ddrestart_cvb,
     >  c,dum,sxc,.false.,cvb,res,rhs,
     >  dum,rhsp,solp,solp_res,
     >  .false.,.false.,.true.,maxdav,nvb,nvb,
     >  nvguess,nvrestart,isaddle,ifollow,mxiter,
     >  resthr,orththr,nortiter,zero,
     >  ioptc,iter,fx,ip(3))
      have_solved_it=.true.

      ovraa=one
      svb=fx
      return
      end
