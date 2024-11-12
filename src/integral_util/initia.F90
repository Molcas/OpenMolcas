!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1995, Roland Lindh                                     *
!***********************************************************************

subroutine InitIA(I,mDeg)
! Purpose: Express the interaction tensor, defined by the
! quantities T(a,b,c) as functions of the vector R=(x,y,z),
! where a,b, and c are nonnegative integers and
! T(a,b,c)=((d/dx)**a)((d/dy)**b)((d/dz)**c) 1/R, in terms
! of a polynomial:
! T(a,b,c)=
!  (sum over p,q,r of I(a,b,c,p,q,r) x**p y**q z**r)/(R**(2*n+1)),
! where n=a+b+c.
! The polynomial coefficients are integers, and are 0 unless
! p+q+r=n.
! Author: PAM

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: mDeg
integer(kind=iwp), intent(out) :: I(0:mDeg,0:mDeg,0:mDeg,0:mDeg,0:mDeg,0:mDeg)
integer(kind=iwp) :: a, b, c, n, new, p, q, r

! initialize:
I(:,:,:,:,:,:) = 0
I(0,0,0,0,0,0) = 1
if (mDeg > 0) then
  I(1,0,0,1,0,0) = -1
  I(0,1,0,0,1,0) = -1
  I(0,0,1,0,0,1) = -1
end if
do n=2,mDeg
  do a=0,n
    do b=0,n-a
      c = n-a-b
      do p=0,n
        do q=0,n-p
          r = n-p-q
          new = 0
          if (a > 0) then
            if (p > 0) new = (p-(2*n))*I(a-1,b,c,p-1,q,r)
            if (q > 1) new = new+(p+1)*I(a-1,b,c,p+1,q-2,r)
            if (r > 1) new = new+(p+1)*I(a-1,b,c,p+1,q,r-2)
          else if (b > 0) then
            if (q > 0) new = (q-(2*n))*I(a,b-1,c,p,q-1,r)
            if (r > 1) new = new+(q+1)*I(a,b-1,c,p,q+1,r-2)
            if (p > 1) new = new+(q+1)*I(a,b-1,c,p-2,q+1,r)
          else
            if (r > 0) new = (r-(2*n))*I(a,b,c-1,p,q,r-1)
            if (p > 1) new = new+(r+1)*I(a,b,c-1,p-2,q,r+1)
            if (q > 1) new = new+(r+1)*I(a,b,c-1,p,q-2,r+1)
          end if
          I(a,b,c,p,q,r) = new
        end do
      end do
    end do
  end do
end do

! write out only elements with a>=b>=c. The others are obtained
! by index permutation.
! This restriction has been removed! (Roland Lindh)
!n = mDeg
!do a=n,0,-1
!  do b=n-a,0,-1
!    c = n-a-b
!    write(u6,'(5x,''T('',i1,'','',i1,'','',i1,'')='',i5)')a,b,c,C_Ind(n,a,c)
!    do p=n,0,-1
!      do q=n-p,0,-1
!        r = n-p-q
!        coef = I(a,b,c,p,q,r)
!        if (coef /= 0) write(u6,'(10x,i8,''*x**'',i1,'' *y**'',i1,'' *z**'',i1,i5)') coef,p,q,r,C_Ind(n,p,r)
!      end do
!    end do
!  end do
!end do

return

end subroutine InitIA
