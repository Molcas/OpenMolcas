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
      subroutine o7b2_cvb(nparm,dx,
     >  dxnrm,grdnrm,close2conv)
      implicit real*8 (a-h,o-z)
      logical skip
      logical close2conv
#include "opt_cvb.fh"
#include "locopt1_cvb.fh"
#include "locopt2_cvb.fh"
#include "trst_cvb.fh"
#include "tune_cvb.fh"

#include "opt2_cvb.fh"
      external asonc7_cvb,ddres2upd7_cvb
      dimension dx(nparm)
      save one,resthr_old
      data one/1d0/

      if(.not.close2conv)then
        resthr_use=1d-5
      else
        resthr_use=5d-2*grdnrm
        resthr_use=min(1d-5,resthr_use)
        resthr_use=max(1d-9,resthr_use)
      endif
      skip=(resthr_use.eq.resthr_old.and.have_solved_it)
      resthr_old=resthr_use
      if(skip)goto 100
      call axex_cvb(asonc7_cvb,ddres2upd7_cvb,dx,
     >  resthr_use,ioptc2,iter2,fx_exp)
      have_solved_it=.true.
      exp=.5d0*fx_exp

      if(ip.ge.2)write(6,'(2a,i4)')' Number of iterations for ',
     >  'direct diagonalization :',iter2

      if(ioptc2.ne.0)then
        write(6,*)' Direct diagonalization not converged!'
        call abend_cvb()
      endif

      if(ip.ge.2)then
        write(6,'(a)')' Eigenvector to be followed :'
        call vecprint_cvb(dx,nparm+1)
      endif
      if(abs(dx(1)).gt.1d-8)then
        fac1=one/dx(1)
      else
        fac1=sign(one,dx(1))
      endif
      call dscal_(nparm,fac1,dx(1),1)
      do 200 i=1,nparm
      dx(i)=dx(i+1)
200   continue
100   dxnrm=dnrm2_(nparm,dx,1)
      if(.not.close2conv)then
        ipu=1
      else
        ipu=2
      endif
      if(dxnrm.gt.hh.or.scalesmall(ipu))then
        call dscal_(nparm,hh/dxnrm,dx,1)
        dxnrm=hh
      endif
      return
      end
