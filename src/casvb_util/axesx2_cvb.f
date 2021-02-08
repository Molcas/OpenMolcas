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
      subroutine axesx2_cvb(asonc,ddres2upd,vec,
     >  resthr_inp,ioptc,iter,fx_exp,
     >  c,axc,sxc,share,res,
     >  ap,solp,solp_res)
      implicit real*8 (a-h,o-z)
      logical share
#include "direct_cvb.fh"
      external asonc,ddres2upd
      external axexsol_cvb,axesxres_cvb,ddrestart_cvb
      dimension c(nparm,maxd),axc(nparm,maxd),sxc(nparm,maxd)
      dimension vec(nparm),res(nparm)
      dimension ap(maxd,maxd)
      dimension solp(maxd),solp_res(maxd),dum(1)

c  If RESTHR_INP unset use default:
      if(resthr_inp.ne.0d0)then
        resthr_use=resthr_inp
      else
        resthr_use=resthr
      endif
      call dirdiag_cvb(asonc,axexsol_cvb,
     >  axesxres_cvb,ddres2upd,ddrestart_cvb,
     >  c,axc,sxc,share,vec,res,dum,
     >  ap,dum,solp,solp_res,
     >  .false.,.true.,.false.,maxd,nparm,nfrdim,
     >  nvguess,nvrestart,isaddle,ifollow,mxit,
     >  resthr_use,orththr,nortiter,corenrg,
     >  ioptc,iter,fx_exp,ip)
      return
      end
