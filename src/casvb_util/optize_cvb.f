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
      subroutine optize_cvb(fx,ioptc,iter,
     >  imetinp,isadinp,mxiter,maxinp,corenrg,ipinp,ipdd1,ipdd2,
     >  strucopt)
c  *********************************************************************
c  *                                                                   *
c  *  Routine for second-order optimization.                           *
c  *                                                                   *
c  *  Uses the following functions:                                    *
c  *                                                                   *
c  *   GETFREE:   Determine # free parameters                          *
c  *    UPDATE:   Add dx to vector of free variables, x                *
c  *        FX:   Calculate f(x)                                       *
c  *      GRAD:   Calculate gradient of f                              *
c  *      HESS:   Calculate the action of the hessian on trial vector  *
c  *   GETHESS:   Calculate full hessian                               *
c  *                                                                   *
c  *  IOPTC is optimization control :                                  *
c  *                                                                   *
c  *  IOPTC=-3    Opt. terminated close to convergence (at request)    *
c  *  IOPTC=-2    Optimization failed -- too small step size           *
c  *  IOPTC=-1    Maximum number of iterations used                    *
c  *  IOPTC= 0    Converged                                            *
c  *  IOPTC= 1    Not complete                                         *
c  *                                                                   *
c  *********************************************************************
      implicit real*8 (a-h,o-z)
      logical maxinp,iter_is_1,strucopt
#include "malloc_cvb.fh"
#include "opt_cvb.fh"

#include "locopt1_cvb.fh"
#include "locopt2_cvb.fh"
#include "trst_cvb.fh"
#include "tune_cvb.fh"

#include "opt2_cvb.fh"
      external dum_a_cvb,o123a_cvb,o123b_cvb,o5b_cvb,o7a_cvb,o7b_cvb
      external o8b_cvb,o10a_cvb,o10b_cvb
      external o12sa_cvb,o12sb_cvb,o12ea_cvb,o12eb_cvb
      save zero
      data zero/0.d0/

      if(mxiter.eq.0)then
        ioptc=-1
        return
      endif

c  Initialize for new optimization - input parameters:
      imethod=imetinp
      isaddle=isadinp
      maxize=maxinp
      ip=ipinp

c  Parameters initialized:
      exp=zero
      hh=hhstart
      hhkeep=hh
      ioptc=1

      if(maxize)then
        ifollow=1
      else
        ifollow=2
      endif

      call fx_cvb(fx,.false.)
      fxbest=fx
      do 1000 iter=1,mxiter
      iter_is_1=(iter.eq.1)
      call getfree_cvb(nparm,n_div,nfrdim,iter,fx)
      if(nfrdim.le.0)then
        if(ip.ge.0)write(6,'(a)')' No free parameters!'
        if(ip.ge.0)write(6,'(a)')' Calculation converged.'
        ioptc=0
        return
      endif
      if(fx.lt.zero.and.maxize)then
        call chgsgn_cvb(fx)
        call getfree_cvb(nparm,n_div,nfrdim,iter,fx)
      endif

      if(imethod.eq.1.or.imethod.eq.2.or.imethod.eq.3)then
        ix(1) = mstackr_cvb(nparm)
        ix(2) = mstackr_cvb(nparm)
        ix(3) = mstackr_cvb(nparm*nparm)
        ix(4) = mstackr_cvb(nparm)
        ix(5) = mstackr_cvb(nparm)
        ix(6) = mstackr_cvb(nparm)
        ix(7) = mstackr_cvb(nparm)
        call optize2_cvb(fx,nparm,ioptc,
     >    w(ix(1)),w(ix(2)),iter_is_1,
     >    o123a_cvb,o123b_cvb)
        call mfreer_cvb(ix(1))
      elseif(imethod.eq.5)then
        ix(1) = mstackr_cvb(nparm)
        ix(2) = mstackr_cvb(nparm)
        call optize2_cvb(fx,nparm,ioptc,
     >    w(ix(1)),w(ix(2)),iter_is_1,
     >    dum_a_cvb,o5b_cvb)
        call mfreer_cvb(ix(1))
      elseif(imethod.eq.7)then
        ix(1) = mstackr_cvb(nparm+1)
        ix(2) = mstackr_cvb(nparm+1)
        maxd=min(nparm+1,200)
        mxit=500
        call ddinit_cvb('AxEx',nparm+1,nfrdim+1,maxd,mxit,
     >    ifollow,isaddle,ipdd1,zero,n_div)
        call asonC7init_cvb(ix(2),ipdd2)
        call optize2_cvb(fx,nparm,ioptc,
     >    w(ix(1)),w(ix(2)),iter_is_1,
     >    o7a_cvb,o7b_cvb)
        call mfreer_cvb(ix(1))
      elseif(imethod.eq.8)then
        ix(1) = mstackr_cvb(nparm)
        ix(2) = mstackr_cvb(nparm)
        ix(3) = mstackr_cvb((nparm+1)*(nparm+1))
        ix(4) = mstackr_cvb(nparm+1)
        call optize2_cvb(fx,nparm,ioptc,
     >    w(ix(1)),w(ix(2)),iter_is_1,
     >    dum_a_cvb,o8b_cvb)
        call mfreer_cvb(ix(1))
      elseif(imethod.eq.9)then
        i1 = mstackr_cvb(nparm)
        i2 = mstackr_cvb(nparm)
        i3 = mstackr_cvb(nparm)
        call optize9_cvb(fx,nparm,ioptc,
     >    w(i1),w(i2),w(i3))
      call mfreer_cvb(i1)
      elseif(imethod.eq.10)then
        ix(1)  = mstackr_cvb(nparm)
        ix(2)  = mstackr_cvb(nparm)
        maxd=min(nparm,200)
        mxit=500
        call ddinit_cvb('AxExb',nparm,nfrdim,maxd,mxit,
     >    ifollow,isaddle,ipdd1,zero,n_div)
        call asonc10init_cvb(ipdd2)
        call optize2_cvb(fx,nparm,ioptc,
     >    w(ix(1)),w(ix(2)),iter_is_1,
     >    o10a_cvb,o10b_cvb)
        call mfreer_cvb(ix(1))
      elseif(imethod.eq.12.and.maxize)then
        if(strucopt)then
          nparm_dav=nparm
          nfrdim_dav=nfrdim
        else
          nparm_dav=nparm+1
          nfrdim_dav=nfrdim+1
        endif
        ix(1)  = mstackr_cvb(nparm_dav)
        ix(2)  = mstackr_cvb(nparm_dav)
        maxd=min(nparm_dav,200)
        mxit=500
        call ddinit_cvb('Axb',nparm_dav,nfrdim_dav,maxd,mxit,
     >    ifollow,isaddle,ipdd1,zero,0)
        call asonc12sinit_cvb(ipdd2)
        call optize2_cvb(fx,nparm_dav,ioptc,
     >    w(ix(1)),w(ix(2)),iter_is_1,
     >    o12sa_cvb,o12sb_cvb)
        call mfreer_cvb(ix(1))
      elseif(imethod.eq.12.and..not.maxize)then
        if(strucopt)then
          nparm_dav=nparm
          nfrdim_dav=nfrdim
        else
          nparm_dav=nparm+1
          nfrdim_dav=nfrdim+1
        endif
        ix(1) = mstackr_cvb(nparm_dav)
        ix(2) = mstackr_cvb(nparm_dav)
        maxd=min(nparm_dav,200)
        mxit=500
        call ddinit_cvb('AxESx',nparm_dav,nfrdim_dav,maxd,mxit,
     >    ifollow,isaddle,ipdd1,corenrg,n_div)
        call asonc12einit_cvb(ipdd2)
        call optize2_cvb(fx,nparm_dav,ioptc,
     >    w(ix(1)),w(ix(2)),iter_is_1,
     >    o12ea_cvb,o12eb_cvb)
        call mfreer_cvb(ix(1))
      else
        write(6,*)' Unrecognized optimization algorithm!',imethod
        call abend_cvb()
      endif
      if(ioptc.le.0)goto 1100
1000  continue
      iter=iter-1
1100  if(ioptc.ge.0)call getfree_cvb(nparm,n_div,nfrdim,-1,fx)
      if(iter.eq.mxiter.and.ioptc.gt.0)ioptc=-1
      return
      end
