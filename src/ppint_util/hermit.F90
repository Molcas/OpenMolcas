!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      subroutine hermit_molcas(nn,x,a,eps)
!                  calculates the zeros  x(i)  of the nn-th order
!                hermite polynomial.  the largest zero will be
!                stored in x(1).  also calculates the corresponding
!                coefficients  a(i)  of the nn-th order gauss-hermite
!                quadrature formula of degree 2*nn-1.  the factor of
!                sqrt(pi) has been removed from the a(i).
!  a. h. stroud & d. secrest, gaussian quadrature formulas,
!  prentice-hall, 1966
!
      implicit real*8 (a-h,o-z)
      parameter (sixth = 1.0d0/6.0d0)
      dimension x(*), a(*)
!
      fn = DBLE(nn)
      n1 = nn - 1
      n2 = (nn+1)/2
      cc = 1.0d0
      s  = 0.0d0
      do 15 i=1,n1
        s  = s + 0.5d0
        cc = s*cc
   15 continue
      s  = (2.0d0*fn+1.0d0)**sixth
      do 30 i=1,n2
        if( i.eq.1 ) then
!         # largest zero
          xt = s**3 - 1.85575d0/s
        elseif( i.eq.2 ) then
!         # second zero
          xt = xt - 1.14d0*fn**0.426d0/xt
        elseif( i.eq.3 ) then
!         # third zero
          xt = 1.86d0*xt - 0.86d0*x(1)
        elseif( i.eq.4 ) then
!         # fourth zero
          xt = 1.91d0*xt - 0.91d0*x(2)
        else
!         # all other zeros
          xt = 2.0d0*xt - x(i-2)
        endif
!
        call hroot(xt,nn,dpn,pn1,eps)
        x(i) = xt
        a(i) = cc/dpn/pn1
!       write (6,'(2i4,2d25.17)') nn, i, xt, a(i)
        ni = nn-i+1
        x(ni) = -xt
        a(ni) = a(i)
   30 continue
      return
      end
