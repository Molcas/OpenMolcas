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
! This routine diagonalize a tridiagonal symmetric matrix using the    *
! QL algorithm. The matrix is stored in lower triangular form.         *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund university, Sweden                                     *
! Written September 2005                                               *
!                                                                      *
!***********************************************************************

implicit none
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
integer n
integer nv
integer irc
real*8 H(*)
real*8 U(nv,n)
!----------------------------------------------------------------------*
! Parameters                                                           *
! MxDim - Largest case that can be handled                             *
! zThr  - Threshold for when an element is regarded as zero            *
!----------------------------------------------------------------------*
integer MxDim
parameter(MxDim=5000)
real*8 zThr
parameter(zThr=1.0d-16)
real*8 qThr
parameter(qThr=1.0d-20)
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
real*8 d(MxDim), e(MxDim)
real*8 g, r, c, s, p, f, b
integer i, j, k, l, m
integer iter, maxiter

!----------------------------------------------------------------------*
! Setup                                                                *
!----------------------------------------------------------------------*
irc = 0
if (n >= MxDim) then
  irc = 1
  !write(6,*) 'QLdiag: system too large!'
  return
end if
if (n < 1) then
  write(6,*) 'QLdiag: zero size system!'
  call Abend()
end if
!----------------------------------------------------------------------*
! Make local copies of diagonal and off-diagonal                       *
!----------------------------------------------------------------------*
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
e(n) = 0.0d0
!----------------------------------------------------------------------*
! Solve it                                                             *
!----------------------------------------------------------------------*
maxiter = 0
do l=1,n
  iter = 0
1 continue
  do j=l,n-1
    if (abs(e(j)) < zThr) then
      m = j
      goto 2
    end if
  end do
  m = n
2 continue
  if (m /= l) then
    if (iter == 25) then
      irc = 1
      !write(6,*) 'QLdiag: ran out of iterations'
      goto 900
    end if
    iter = iter+1
    maxiter = max(maxiter,iter)
    g = (d(l+1)-d(l))/(2.0d0*e(l))
    r = sqrt(1.0d0+g*g)
    g = d(m)-d(l)+e(l)/(g+sign(r,g))
    s = 1.0d0
    c = 1.0d0
    p = 0.0d0
    do i=m-1,l,-1
      f = s*e(i)
      b = c*e(i)
      r = sqrt(f*f+g*g)
      e(i+1) = r
      if (abs(r) <= qThr) then
        d(i+1) = d(i+1)-p
        e(m) = 0.0d0
        goto 1
      end if
      s = f/r
      c = g/r
      g = d(i+1)-p
      r = (d(i)-g)*s+2.0d0*c*b
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
    e(m) = 0.0d0
    goto 1
  end if
end do
!----------------------------------------------------------------------*
! Copy back local copy                                                 *
!----------------------------------------------------------------------*
900 continue
!write(6,*) 'QLdiag: maxiter ',maxiter
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

!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine QLdiag
