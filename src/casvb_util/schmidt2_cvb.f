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
      subroutine schmidt2_cvb(c,sxc,cnrm,nvec,sao,n,metr)
      implicit real*8 (a-h,o-z)
      dimension c(n,nvec),sxc(n,nvec),cnrm(nvec),sao(*)
      save thresh
      data thresh/1d-20/

      do 100 i=1,nvec
      do 200 j=1,i-1
      if(cnrm(j).gt.thresh)
     >  call daxpy_(n,-ddot_(n,c(1,i),1,sxc(1,j),1)/cnrm(j),
     >  c(1,j),1,c(1,i),1)
200   continue
      if(metr.ne.0)call saoon_cvb(c(1,i),sxc(1,i),1,sao,n,metr)
      cnrm(i)=ddot_(n,c(1,i),1,sxc(1,i),1)
100   continue
      return
      end
