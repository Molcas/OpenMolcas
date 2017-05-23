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
      subroutine mxsqrt_cvb(a,n,ipow)
      implicit real*8 (a-h,o-z)
#include "malloc_cvb.fh"
      dimension a(n,n)

      i1 = mstackr_cvb(n)
      i2 = mstackr_cvb(n*n)
      i3 = mstackr_cvb(n)
      i4 = mstackr_cvb(n)
      i5 = mstackr_cvb(n*n)
      ifail=0
      call casvb_rs(n,n,a,w(i1),1,w(i2),w(i3),w(i4),ifail)
      if(ifail.ne.0)then
        write(6,*)' Fatal error in diagonalization (MXSQRT) :',ifail
        call abend_cvb()
      endif
      call fzero(a,n*n)
      do 100 i=1,n
100   a(i,i)=sqrt(w(i+i1-1))**ipow
      call mxatb_cvb(w(i2),a,n,n,n,w(i5))
      call fzero(a,n*n)
      do 200 k=1,n
      do 200 j=1,n
      do 200 i=1,n
200   a(i,j)=a(i,j)+w(i+(k-1)*n+i5-1)*w(j+(k-1)*n+i2-1)
      call mfreer_cvb(i1)
      return
      end
