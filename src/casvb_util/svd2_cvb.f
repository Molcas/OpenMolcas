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
      subroutine svd2_cvb(ainp,val,vec,vmat,n1,n2,n12,
     >  a,w,u,v,rv1,indx)
      implicit real*8 (a-h,o-z)
      dimension ainp(n1,n2),val(n2),vec(n1,n2),vmat(n2,n2)
      dimension a(n12,n2),w(n2),u(n12,n2),v(n12,n2),rv1(n2),indx(n2)

      if(n12.eq.n1)then
        call fmove_cvb(ainp,a,n1*n2)
      else
        call fzero(a,n12*n2)
        do 100 i=1,n2
        call fmove_cvb(ainp(1,i),a(1,i),n1)
100     continue
      endif
      ierr=0
      call svd(n12,n1,n2,a,w,.true.,u,.true.,v,ierr,rv1)

      if(ierr.ne.0)then
        write(6,*)' Fatal error in SVD_CVB!',ierr
        call abend_cvb()
      endif

c  Eispack code is broken, in the following u is generated
c  from v :

c  First recreate a :
      if(n12.eq.n1)then
        call fmove_cvb(ainp,a,n1*n2)
      else
        call fzero(a,n12*n2)
        do 200 i=1,n2
        call fmove_cvb(ainp(1,i),a(1,i),n1)
200     continue
      endif

      do 300 i=1,n2
      call mxatb_cvb(a,v(1,i),n12,n2,1,u(1,i))
      call dscal_(n12,1d0/dnrm2_(n12,u(1,i),1),u(1,i),1)
300   continue

c  Sort singular values in ascending order:
      call sortindxr_cvb(n2,w,indx)
      do 400 i=1,n2
      val(i)=w(indx(i))
      call fmove_cvb(v(1,indx(i)),vmat(1,i),n2)
      call fmove_cvb(u(1,indx(i)),vec(1,i),n1)
400   continue
      return
      end
      function detm_cvb(a,n)
      implicit real*8 (a-h,o-z)
#include "malloc_cvb.fh"
      dimension a(n*n)
      dimension det(2)
cstart linpack_determinant
      save zero,one
      data zero/0d0/,one/1d0/
celse
c;      save one
c;      data one/1d0/
cend

      if(n.eq.0)then
        detm_cvb=one
        return
      endif
      i1 = mstackr_cvb(n*n)
      i2 = mstacki_cvb(n)
      ierr=0
      call fmove_cvb(a,w(i1),n*n)
      call dgetrf_(n,n,w(i1),n,iw(i2),ierr)
cstart linpack_determinant
c      call dgefa(w(i1),n,n,iw(i2),ierr)
      i3 = mstackr_cvb(n*n)
      if(ierr.ne.0)then
        detm_cvb=zero
        call mfreer_cvb(i1)
        return
      endif
      call dgedi(w(i1),n,n,iw(i2),det,w(i3),10)
celse
c;      dl=0d0
c;      ds=1d0
c;      do k=0,n-1
c;        dl=dl+log10(abs(w(i1+k*(n+1))))
c;        if(w(i1+k*(n+1)).lt.0d0)ds=-ds
c;      end do
c;      det(2)=dble(int(dl))
c;      det(1)=ds*(10d0**(dl-det(2)))
cend
      detm_cvb=det(1) * 10d0**det(2)
      call mfreer_cvb(i1)
      return
      end
c  *******************************************
c  ** Orthogonalisation, normalisation etc. **
c  *******************************************
