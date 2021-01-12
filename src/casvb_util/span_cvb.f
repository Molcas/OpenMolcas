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
      subroutine span_cvb(a,nvec,nlin,s,n,metr)
c  Creates span of a vector set.
c  On entry : A is NVEC by N
c  On exit  : NLIN is the number of linearly independent vectors
c               (NVEC and NLIN may coincide in calling sequence)
c             A contains set of NLIN linearly independent vectors
c               spanning the same set as the NVEC input vectors
c               Vectors will be orthonormal on exit
      implicit real*8 (a-h,o-z)
      dimension a(n,nvec),s(*)
      save thresh
      data thresh/1.d-10/
      nvect=nvec
      ierr=1
      call nize_cvb(a,nvect,s,n,metr,ierr)
      call schmidt_cvb(a,nvect,s,n,metr)
      nlin=0
      do 100 i=1,nvect
      cnrm=dnrm2_(n,a(1,i),1)
      if(cnrm.gt.thresh)then
        nlin=nlin+1
        call fmove_cvb(a(1,i),a(1,nlin),n)
      endif
100   continue
      ierr=1
      call nize_cvb(a,nlin,s,n,metr,ierr)
      return
      end
