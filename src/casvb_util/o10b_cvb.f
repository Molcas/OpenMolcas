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
      subroutine o10b_cvb(nparm,dxnrm,grdnrm,close2conv)
      implicit real*8 (a-h,o-z)
      logical close2conv
#include "malloc_cvb.fh"
#include "opt_cvb.fh"
#include "locopt1_cvb.fh"
#include "locopt2_cvb.fh"
#include "trst_cvb.fh"
#include "tune_cvb.fh"

#include "opt2_cvb.fh"
      external asonc10_cvb,ddres2upd10_cvb

      if(.not.close2conv)then
        resthr_use=1d-5
      else
        resthr_use=5d-2*grdnrm
        resthr_use=min(1d-5,resthr_use)
        resthr_use=max(1d-9,resthr_use)
      endif
      call axexb_cvb(asonc10_cvb,ddres2upd10_cvb,w(ix(1)),
     >  resthr_use,ioptc2,iter2,fx_exp)
      have_solved_it=.true.

      if(ip.ge.2)write(6,'(2a,i4)')' Number of iterations for ',
     >  'direct diagonalization :',iter2

      if(ioptc2.ne.0)then
        write(6,*)' Direct diagonalization not converged!'
        call abend_cvb()
      endif

      dxnrm=dnrm2_(nparm,w(ix(1)),1)
      return
      end
