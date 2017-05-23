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
      subroutine mxunit_cvb(a,n)
      implicit real*8 (a-h,o-z)
      dimension a(n,n)
      save one
      data one/1d0/

      call fzero(a,n*n)
      do 100 i=1,n
100   a(i,i)=one
      return
      end
      logical function mxorth_cvb(a,n)
c  Returns .TRUE. if A is orthogonal.
      implicit real*8 (a-h,o-z)
#include "malloc_cvb.fh"
      dimension a(n,n)
      save thresh,one
      data thresh/1d-8/,one/1d0/

      i1 = mstackr_cvb(n*n)
      i2 = mstackr_cvb(n*n)
c  W(I1) <= A transpose
      do 100 i=1,n
      do 100 j=1,n
100   w(i+(j-1)*n+i1-1)=a(j,i)
      call mxatb_cvb(w(i1),a,n,n,n,w(i2))
c  W(I2) identity ??
      mxorth_cvb=.true.
      do 200 j=1,n
      do 200 i=1,n
      if(i.ne.j)then
        tst=abs(w(i+(j-1)*n+i2-1))
        if(tst.gt.thresh)mxorth_cvb=.false.
      else
        tst=abs(w(i+(j-1)*n+i2-1)-one)
        if(tst.gt.thresh)mxorth_cvb=.false.
      endif
200   continue
      return
      end
