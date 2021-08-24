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
      subroutine axb_cvb(asonc,ddres2upd,vec,
     >  resthr_inp,ioptc,iter,fx_exp)
c  *********************************************************************
c  *                                                                   *
c  *  DIRDIAG front-end for solving  A x = b .                         *
c  *                                                                   *
c  *********************************************************************
      implicit real*8 (a-h,o-z)
#include "WrkSpc.fh"
#include "direct_cvb.fh"
      external asonc,ddres2upd

      dimension vec(*)

      call axb2_cvb(asonc,ddres2upd,vec,
     >  resthr_inp,ioptc,iter,fx_exp,
     >  work(idd(1)),work(idd(2)),work(idd(3)),work(idd(4)),
     >  work(idd(5)),work(idd(6)),work(idd(7)))
      return
      end
