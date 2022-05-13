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

subroutine Lobatto(ndeg,trw)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Three
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ndeg
real(kind=wp), intent(out) :: trw(3*(ndeg+2)*(ndeg+3)/2)
integer(kind=iwp) :: ii, ir, jj, jr, k, n
real(kind=wp) :: c, delta, dmax, f, fnew, fold, fp, fpnew, fpold, rk, x, xterm
real(kind=wp), allocatable :: recurs(:), roots(:,:), wghts(:,:)

call mma_allocate(roots,ndeg,ndeg,label='roots')
call mma_allocate(recurs,ndeg,label='recurs')
! Start: Accurately known are p0=1 and p1=x
roots(1,1) = Zero
! Recursion coefficients:
do k=1,ndeg
  rk = real(k,kind=wp)
  recurs(k) = (rk*(rk+Two))/((Two*rk+One)*(Two*rk+Three))
end do

do k=2,ndeg
  ! Construct a start approximation to the roots:
  roots(1,k) = (real(k,kind=wp)*(roots(1,k-1)+One))/real(k+1,kind=wp)-One
  roots(k,k) = (real(k,kind=wp)*(roots(k-1,k-1)-One))/real(k+1,kind=wp)+One
  do ir=2,k-1
    roots(ir,k) = (real(k+1-ir,kind=wp)*roots(ir,k-1)+real(ir,kind=wp)*roots(ir-1,k-1))/real(k+1,kind=wp)
  end do
  ! Start modified Newton-Raphson iterations. Parallell treatment of roots:
  do
    dmax = Zero
    do ir=1,k
      ! Compute value and derivative of polynomial:
      x = roots(ir,k)
      fpold = Zero
      fold = One
      fp = One
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
      xterm = Zero
      do jr=1,k
        if (jr /= ir) xterm = xterm+One/(x-roots(jr,k))
      end do
      ! Update:
      delta = -f/(fp-f*xterm)
      roots(ir,k) = roots(ir,k)+delta
      dmax = max(dmax,abs(delta))
    end do
    if (dmax <= 1.0e-12_wp) exit
  end do
end do

call mma_deallocate(recurs)
call mma_allocate(wghts,ndeg,ndeg,label='wghts')

! Compute weights:
do k=1,ndeg
  do ir=1,k
    x = roots(ir,k)
    fold = One
    f = x
    do n=1,k
      fnew = x*f*(Two*real(n,kind=wp)+One)/(real(n,kind=wp)+One)-fold*real(n,kind=wp)/(real(n,kind=wp)+One)
      fold = f
      f = fnew
    end do
    wghts(ir,k) = Two/(f*f*real(k+1,kind=wp)*real(k+2,kind=wp))
  end do
end do

do n=3,ndeg+2
  trw(3*n*(n-1)/2+1) = -One                       ! (n-1,1,1)   n=nDeg+2
  trw(3*n*(n-1)/2+2) = Two/real(n*(n-1),kind=wp)  ! (n-1,1,2)
  trw(3*n*(n+1)/2-2) = One                        ! (n-1,n-1,1)
  trw(3*n*(n+1)/2-1) = Two/real(n*(n-1),kind=wp)  ! (n-1,n-1,2)
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

trw(1:9) = Zero

do k=1,ndeg
  n = k+1
  do ir=1,k
    ii = 3*n*(n+1)/2+ir*3+1  ! (n+1,ir,1)   n=nDeg+1, ir=2,nDeg
    jj = 3*n*(n+1)/2+ir*3+2  ! (n+1,ir,2)
    trw(ii) = roots(ir,k)
    trw(jj) = wghts(ir,k)
  end do
end do

call mma_deallocate(roots)
call mma_deallocate(wghts)

!write(u6,*) 'Lobatto'
!do i=1,ndeg+2
!  write(u6,*) 'i=',i
!  do j=1,i
!    write(u6,*) trw(3*i*(i-1)/2+3*(j-1)+1),trw(3*i*(i-1)/2+3*(j-1)+2),trw(3*i*(i-1)/2+3*(j-1)+3)
!  end do
!end do

return

end subroutine Lobatto
