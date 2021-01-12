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
      subroutine nize_cvb(c,nnrm,s,n,metr,ierr)
c  Normalizes NNRM vectors in C.
      implicit real*8 (a-h,o-z)
      logical safe
#include "malloc_cvb.fh"
      dimension c(n,nnrm),s(*)
      save one,thresh
      data one/1.d0/,thresh/1.d-8/

      if(metr.ne.0)i1=mstackr_cvb(n)
      safe=ierr.ne.0
      do 100 i=1,nnrm
      if(metr.eq.0)then
        cnrm=dnrm2_(n,c(1,i),1)
      else
        call saoon_cvb(c(1,i),w(i1),1,s,n,metr)
        cnrm=sqrt(ddot_(n,c(1,i),1,w(i1),1))
      endif
      if(safe.and.cnrm.lt.thresh)then
        ierr=ierr+1
      else
        call dscal_(n,one/cnrm,c(1,i),1)
      endif
100   continue
      if(metr.ne.0)call mfreer_cvb(i1)
      return
      end
