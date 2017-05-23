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
      subroutine testconv_cvb(fx,npr,
     >  dx,w2,exp_tc,
     >  close2conv,converged,wrongstat)
      implicit real*8 (a-h,o-z)
      logical close2conv,converged,wrongstat
#include "opt_cvb.fh"
#include "locopt1_cvb.fh"
#include "locopt2_cvb.fh"
#include "trst_cvb.fh"
#include "tune_cvb.fh"

      dimension dx(npr),w2(npr)
      save one
      data one/1.d0/

      if(maxize)then
        nposeig=min(isaddle,npr)
      else
        nposeig=max(npr-isaddle,npr)
      endif
      nnegeig=npr-nposeig

c  EIGMX is maximum of NNEGEIG first Hessian eigenvalues (which should
c  all be negative) EIGMN the minimum of NPOSEIG last eigenvalues :
      eigmx=-one
      eigmn=one
      eigmna=one

      call zz_cvb(act,zz,fx,fxbest,exp_tc,ip)
      fxbest=fx

c  << Test for convergence or near-convergence >>
      call testconv2_cvb(close2conv,converged,wrongstat,
     >  act,zz,dx,w2,npr,eigmn,eigmx,eigmna,nposeig,nnegeig)
      return
      end
