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
      subroutine mxinv_cvb(a,n)
      implicit real*8 (a-h,o-z)
#include "malloc_cvb.fh"
      dimension a(n,n)
      save one,thresh
      data one/1.d0/,thresh/1.d-10/

      i1 = mstackr_cvb(n*n)
      i2 = mstackr_cvb(n*n)
      i3 = mstacki_cvb(n)
      ierr=0
      call fmove_cvb(a,w(i1),n*n)
      call dgetrf_(n,n,w(i1),n,iw(i3),ierr)
      if(ierr.ne.0) then
        write(6,*)' Error in LU decomposition for inversion:',ierr
        call mxprint_cvb(a,n,n,0)
        call abend_cvb()
      endif
      call dgetri_(n,w(i1),n,iw(i3),w(i2),n*n,ierr)
c  Check solution
      call mxatb_cvb(a,w(i1),n,n,n,w(i2))
      do 100 i=1,n
      w(i+(i-1)*n+i2-1)=w(i+(i-1)*n+i2-1)-one
100   continue
      rms=sqrt(ddot_(n*n,w(i2),1,w(i2),1)/DBLE(n*n))
      if(rms.gt.thresh)then
        write(6,*)' Fatal error in matrix inversion - error:',rms
        write(6,*)' Singular or near-singular matrix.'
        write(6,*)' Matrix :'
        call mxprint_cvb(a,n,n,0)
        write(6,*)' Inverted matrix :'
        call mxprint_cvb(w(i1),n,n,0)
        write(6,*)' Check :'
        call mxprint_cvb(w(i2),n,n,0)
        call mxdiag_cvb(a,w(i2),n)
        write(6,*)' Eigenvalues :'
        call mxprint_cvb(w(i2),1,n,0)
        write(6,*)' Eigenvectors :'
        call mxprint_cvb(a,1,n,0)
        call abend_cvb()
      endif
      call fmove_cvb(w(i1),a,n*n)
      call mfreer_cvb(i1)
      return
      end
