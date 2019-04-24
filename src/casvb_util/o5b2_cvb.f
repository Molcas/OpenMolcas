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
      subroutine o5b2_cvb(nparm,
     >  dx,grad,
     >  dxnrm,close2conv)
      implicit real*8 (a-h,o-z)
      logical close2conv
#include "opt_cvb.fh"
#include "locopt1_cvb.fh"
#include "locopt2_cvb.fh"
#include "trst_cvb.fh"
#include "tune_cvb.fh"

      dimension dx(nparm),grad(nparm)
      save one
      data one/1d0/

      call fmove_cvb(grad,dx,nparm)
      if(.not.maxize)call dscal_(nparm,-one,dx,1)
      dxnrm=dnrm2_(nparm,dx,1)
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
