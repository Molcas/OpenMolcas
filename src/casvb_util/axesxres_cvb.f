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
      subroutine axesxres_cvb(axc,sxc,dum,res,
     >  solp_res,maxdav,n,itdav,eig_res,
     >  is_converged)
      implicit real*8 (a-h,o-z)
      logical is_converged
#include "direct_cvb.fh"
      dimension axc(n,maxdav),sxc(n,maxdav),res(n)
      dimension solp_res(maxdav)

      call fzero(res,n)
      do 300 i=1,itdav
      do 301 ivb=1,n
      res(ivb)=res(ivb)+(axc(ivb,i)-eig_res*sxc(ivb,i))*solp_res(i)
301   continue
300   continue

      is_converged=(jroot.eq.iroot)
      return
c Avoid unused argument warnings
      if (.false.) call Unused_real(dum)
      end
