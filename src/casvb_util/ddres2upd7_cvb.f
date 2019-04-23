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
      subroutine ddres2upd7_cvb(res,c,n)
      implicit real*8 (a-h,o-z)
#include "direct_cvb.fh"
      dimension res(n),c(n)

      if(n_div.eq.0)then
        call fmove_cvb(res,c,n)
      else
        resnrm1=dnrm2_(n_div-1,res(2),1)
        resnrm2=dnrm2_(n-n_div,res(n_div+1),1)
        if(resnrm1.gt.resnrm2)then
          call fmove_cvb(res,c,n_div)
          call fzero(c(n_div+1),n-n_div)
        else
          call fzero(c,n_div)
          c(1)=res(1)
          call fmove_cvb(res(n_div+1),c(n_div+1),n-n_div)
        endif
      endif
      call ddproj_cvb(c,n)
      return
      end
