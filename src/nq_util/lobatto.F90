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

subroutine Lobatto(ndeg,Trw)

implicit real*8(a-h,o-z)
parameter(mxdeg=100)
dimension roots(mxdeg,mxdeg), wghts(mxdeg,mxdeg)
dimension recurs(mxdeg), trw(3*(ndeg+2)*(ndeg+3)/2)

! Start: Accurately known are p0=1 and p1=x
roots(1,1) = 0.0d0
! Recursion coefficients:
do k=1,ndeg
  rk = dble(k)*1.0d0
  recurs(k) = (rk*(rk+2.0d0))/((2.0d0*rk+1.0d0)*(2.0d0*rk+3.0d0))
end do

do k=2,ndeg
  ! Construct a start approximation to the roots:
  roots(1,k) = (dble(k)*(roots(1,k-1)+1.0d0))/dble(k+1)-1.0d0
  roots(k,k) = (dble(k)*(roots(k-1,k-1)-1.0d0))/dble(k+1)+1.0d0
  do ir=2,k-1
    roots(ir,k) = (dble(k+1-ir)*roots(ir,k-1)+dble(ir)*roots(ir-1,k-1))/dble(k+1)
  end do
  ! Start modified Newton-Raphson iterations. Parallell treatment of roots:
  do
    dmax = 0.0d0
    do ir=1,k
      ! Compute value and derivative of polynomial:
      x = roots(ir,k)
      fpold = 0.0d0
      fold = 1.0d0
      fp = 1.0d0
      f = x
      do n=2,k
        c = recurs(n-1)
        fpnew = x*fp+f-c*fpold
        fnew = x*f-c*fold
        fpold = fp
        fold = f
        fp = fpnew
        f = fnew
      end do
      ! Compute the extra denominator term:
      xterm = 0.0d0
      do jr=1,k
        if (jr /= ir) xterm = xterm+1.0d0/(x-roots(jr,k))
      end do
      ! Update:
      delta = -f/(fp-f*xterm)
      roots(ir,k) = roots(ir,k)+delta
      dmax = max(dmax,abs(delta))
    end do
    if (dmax <= 1.0d-12) exit
  end do
end do

! Compute weights:
do k=1,ndeg
  do ir=1,k
    x = roots(ir,k)
    fold = 1.0d0
    f = x
    do n=1,k
      fnew = x*f*(2.0d0*dble(n)+1.0d0)/(dble(n)+1.0d0)-fold*dble(n)/(dble(n)+1.0d0)
      fold = f
      f = fnew
    end do
    wghts(ir,k) = 2.0d0/(f*f*dble(k+1)*dble(k+2))
  end do
end do

do n=3,ndeg+2
  trw(3*n*(n-1)/2+1) = -1.0d0               ! (n-1,1,1)   n=nDeg+2
  trw(3*n*(n-1)/2+2) = 2.0d0/dble(n*(n-1))  ! (n-1,1,2)
  trw(3*n*(n+1)/2-2) = 1.0d0                ! (n-1,n-1,1)
  trw(3*n*(n+1)/2-1) = 2.0d0/dble(n*(n-1))  ! (n-1,n-1,2)
end do

! (1,1,1)
! (1,1,2)
! (1,1,3)
! (2,1,1)
! (2,1,1)
! (2,1,2)
! (2,2,3)
! (2,2,2)
! (2,2,3)

do i=1,9
  trw(i) = 0.0d0
end do

do k=1,ndeg
  n = k+1
  do ir=1,k
    ii = 3*n*(n+1)/2+ir*3+1  ! (n+1,ir,1)   n=nDeg+1, ir=2,nDeg
    jj = 3*n*(n+1)/2+ir*3+2  ! (n+1,ir,2)
    trw(ii) = roots(ir,k)
    trw(jj) = wghts(ir,k)
  end do
end do

!write(6,*) 'Lobatto'
!do i=1,ndeg+2
!  write(6,*) 'i=',i
!  do j=1,i
!    write(6,*) trw(3*i*(i-1)/2+3*(j-1)+1),trw(3*i*(i-1)/2+3*(j-1)+2),trw(3*i*(i-1)/2+3*(j-1)+3)
!  end do
!end do

return

end subroutine Lobatto
