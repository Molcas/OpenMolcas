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
      subroutine axexb_cvb(asonc,ddres2upd,vec,
     >  resthr_inp,ioptc,iter,fx_exp)
c  *********************************************************************
c  *                                                                   *
c  *  DIRDIAG front-end for solving  A x - E x = b .                   *
c  *                                                                   *
c  *********************************************************************
      implicit real*8 (a-h,o-z)
#include "malloc_cvb.fh"
#include "direct_cvb.fh"
      external asonc,ddres2upd

      dimension vec(*)

      call axesxb2_cvb(asonc,ddres2upd,vec,
     >  resthr_inp,ioptc,iter,fx_exp,
     >  w(idd(1)),w(idd(2)),w(idd(1)),.true.,w(idd(3)),w(idd(4)),
     >  w(idd(5)),w(idd(6)),w(idd(7)),w(idd(8)))
      return
      end
