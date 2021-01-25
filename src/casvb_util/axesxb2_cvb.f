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
      subroutine axesxb2_cvb(asonc,ddres2upd,vec,
     >  resthr_inp,ioptc,iter,fx_exp,
     >  c,axc,sxc,share,res,rhs,
     >  ap,rhsp,solp,solp_res)
      implicit real*8 (a-h,o-z)
      logical share
#include "direct_cvb.fh"
      external asonc,ddres2upd
      external axexbsol_cvb,axesxbres_cvb,ddrestart_cvb
      dimension c(nparm,maxd),axc(nparm,maxd),sxc(nparm,maxd)
      dimension vec(nparm),res(nparm),rhs(*)
      dimension ap(maxd,maxd),rhsp(maxd)
      dimension solp(maxd),solp_res(maxd)

c  If RESTHR_INP unset use default:
      if(resthr_inp.ne.0d0)then
        resthr_use=resthr_inp
      else
        resthr_use=resthr
      endif
c------------------------------------------
c      ASonC             ASonC10_CVB
c      AxExbSOL_CVB      DDSOL10_CVB
c      AxESxbRES_CVB     DDRES10_CVB
c      DDRES2UPD         DDRES2UPD10_CVB
c      DDRESTART_CVB     DDRESTART_CVB
c------------------------------------------
      call dirdiag_cvb(asonc,axexbsol_cvb,
     >  axesxbres_cvb,ddres2upd,ddrestart_cvb,
     >  c,axc,sxc,share,vec,res,rhs,
     >  ap,rhsp,solp,solp_res,
     >  .false.,.true.,.true.,maxd,nparm,nfrdim,
     >  nvguess,nvrestart,isaddle,ifollow,mxit,
     >  resthr_use,orththr,nortiter,corenrg,
     >  ioptc,iter,fx_exp,ip)
      return
      end
