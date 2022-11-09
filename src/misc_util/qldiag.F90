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
! Copyright (C) 2005, Per-Olof Widmark                                 *
!***********************************************************************

subroutine QLdiag(H,U,n,nv,irc)
!***********************************************************************
!                                                                      *
! This routine diagonalizes a tridiagonal symmetric matrix using the   *
! QL algorithm. The matrix is stored in lower triangular form.         *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund university, Sweden                                     *
! Written September 2005                                               *
!                                                                      *
!***********************************************************************
!----------------------------------------------------------------------*
! Dummy arguments                                                      *
! n   - Dimension of matrix                                            *
! nv  - Length of eigenvectors nv>=n                                   *
! H   - Matrix to be diagonalized                                      *
! U   - Eigenvectors                                                   *
! irc - Return code, 0 = OK                                            *
!                    1 = Not converged                                 *
!                    2 = Too large system.                             *
!----------------------------------------------------------------------*
! Parameters                                                           *
! MxDim - Largest case that can be handled                             *
! zThr  - Threshold for when an element is regarded as zero            *
!----------------------------------------------------------------------*

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: n, nv
real(kind=wp), intent(inout) :: H(*), U(nv,n)
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp) :: i, iter, j, k, l, m, maxiter
real(kind=wp) :: b, c, f, g, p, r, s
real(kind=wp), allocatable :: d(:), e(:)
integer(kind=iwp), parameter :: MxDim = 5000
real(kind=wp), parameter :: qThr = 1.0e-20_wp, zThr = 1.0e-16_wp

!----------------------------------------------------------------------*
! Setup                                                                *
!----------------------------------------------------------------------*
irc = 0
if (n >= MxDim) then
  irc = 1
  !write(u6,*) 'QLdiag: system too large!'
  return
end if
if (n < 1) then
  write(u6,*) 'QLdiag: zero size system!'
  call Abend()
end if
!----------------------------------------------------------------------*
! Make local copies of diagonal and off-diagonal                       *
!----------------------------------------------------------------------*
call mma_allocate(d,n,label='d')
call mma_allocate(e,n,label='e')
j = 1
do i=1,n
  d(i) = H(j)
  j = j+i+1
end do
j = 2
do i=1,n-1
  e(i) = H(j)
  j = j+i+2
end do
e(n) = Zero
!----------------------------------------------------------------------*
! Solve it                                                             *
!----------------------------------------------------------------------*
maxiter = 0
outer: do l=1,n
  iter = 0
  inner: do
    m = n
    do j=l,n-1
      if (abs(e(j)) < zThr) then
        m = j
        exit
      end if
    end do
    if (m == l) exit inner
    if (iter == 25) then
      irc = 1
      !write(u6,*) 'QLdiag: ran out of iterations'
      exit outer
    end if
    iter = iter+1
    maxiter = max(maxiter,iter)
    g = (d(l+1)-d(l))/(Two*e(l))
    r = sqrt(One+g*g)
    g = d(m)-d(l)+e(l)/(g+sign(r,g))
    s = One
    c = One
    p = Zero
    do i=m-1,l,-1
      f = s*e(i)
      b = c*e(i)
      r = sqrt(f*f+g*g)
      e(i+1) = r
      if (abs(r) <= qThr) then
        d(i+1) = d(i+1)-p
        e(m) = Zero
        cycle inner
      end if
      s = f/r
      c = g/r
      g = d(i+1)-p
      r = (d(i)-g)*s+Two*c*b
      p = s*r
      d(i+1) = g+p
      g = c*r-b
      do k=1,nv
        f = U(k,i+1)
        U(k,i+1) = s*U(k,i)+c*f
        U(k,i) = c*U(k,i)-s*f
      end do
    end do
    d(l) = d(l)-p
    e(l) = g
    e(m) = Zero
  end do inner
end do outer
!----------------------------------------------------------------------*
! Copy back local copy                                                 *
!----------------------------------------------------------------------*
!write(u6,*) 'QLdiag: maxiter ',maxiter
j = 1
do i=1,n
  H(j) = d(i)
  j = j+i+1
end do
j = 2
do i=1,n-1
  H(j) = e(i)
  j = j+i+2
end do
call mma_deallocate(d)
call mma_deallocate(e)

!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine QLdiag
