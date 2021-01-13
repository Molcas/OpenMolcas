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
      subroutine o12eb2_cvb(orbs,cvb,nparm1,nvb,
     >  nfrorb,
     >  gjorb,gjorb2,gjorb3,
     >  dx,
     >  dxnrm,grdnrm,close2conv,strucopt)
      implicit real*8 (a-h,o-z)
      logical strucopt,skip
      logical close2conv
      external asonc12e_cvb,ddres2upd10_cvb
#include "malloc_cvb.fh"
#include "opt_cvb.fh"
#include "locopt1_cvb.fh"
#include "locopt2_cvb.fh"
#include "trst_cvb.fh"
#include "tune_cvb.fh"

#include "davtune_cvb.fh"
#include "opt2_cvb.fh"
      dimension orbs(*),cvb(nvb)
      dimension gjorb(*),gjorb2(*),gjorb3(*)
      dimension dx(nparm1)
      save resthr_old,one
      data one/1d0/

      if(.not.close2conv)then
        resthr_use=1d-5
      else
        resthr_use=5d-2*grdnrm
        resthr_use=max(.5d-5*.6d0,resthr_use)
        resthr_use=min(1d-5,resthr_use)
      endif
      skip=(resthr_use.eq.resthr_old.and.have_solved_it)
      resthr_old=resthr_use
      if(skip)goto 100

      call makegjorbs_cvb(orbs,gjorb,gjorb2,gjorb3)

      call axesx_cvb(asonc12e_cvb,ddres2upd10_cvb,dx,
     >  resthr_use,ioptc2,iter2,fx_exp)
      exp=fx_exp-fxbest
      have_solved_it=.true.

      if(ip.ge.2)write(6,'(2a,i4)')' Number of iterations for ',
     >  'direct diagonalization :',iter2

      if(strucopt)then
        cnrm2=ddot_(nvb,cvb,1,dx(nfrorb+1),1)
c  "Orthogonalize" on CVB to get smallest possible update norm :
        call daxpy_(nvb,-cnrm2,cvb,1,dx(nfrorb+1),1)
c  Scale variables according to overlap with CVB :
        call dscal_(nparm1,one/cnrm2,dx,1)
      else
c  We are doing "Augmented" calc:
        fac=one/dx(1)
c  Scale variables according to overlap with CVB :
        do 50 i=1,nparm1-1
        dx(i)=fac*dx(i+1)
50      continue
      endif

100   dxnrm=dnrm2_(nparm1,dx,1)
      if(.not.close2conv)then
        ipu=1
      else
        ipu=2
      endif
      if(dxnrm.gt.hh.or.scalesmall(ipu))then
        call dscal_(nparm1,hh/dxnrm,dx,1)
        dxnrm=hh
      endif
      return
      end
