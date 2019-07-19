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
      subroutine o10a_cvb(nparm1)
      implicit real*8 (a-h,o-z)
#include "malloc_cvb.fh"
#include "direct_cvb.fh"
#include "opt2_cvb.fh"

      call ddnewopt_cvb()
      have_solved_it=.false.

      ixp=mstackr_cvb(nparm)
      call fmove_cvb(w(ix(2)),w(ixp),nparm)
      call ddproj_cvb(w(ixp),nparm)
      cnrm1=dnrm2_(n_div,w(ixp),1)
      cnrm2=dnrm2_(nparm-n_div,w(n_div+ixp),1)
      if(cnrm1.gt.cnrm2)then
        call ddguess_cvb(w(ixp),n_div,0)
        if(cnrm2.gt.1d-8)
     >    call ddguess_cvb(w(n_div+ixp),nparm-n_div,n_div)
      else
        call ddguess_cvb(w(n_div+ixp),nparm-n_div,n_div)
        if(cnrm1.gt.1d-8)call ddguess_cvb(w(ixp),n_div,0)
      endif
      call ddrhs_cvb(w(ixp),nparm,0)
      call mfreer_cvb(ixp)
      return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(nparm1)
      end
