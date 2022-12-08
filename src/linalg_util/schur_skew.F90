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
! Copyright (C) 2021, Ignacio Fdez. Galvan                             *
!***********************************************************************
!  Schur_Skew
!
!> @brief Compute the real Schur decomposition of an antisymmetric (skew-symmetric) matrix
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Computes the real Schur decomposition of an antisymmetric real matrix \f$ A \f$.
!> The real Schur decomposition satisfies \f$ A = Z E Z^T \f$ with \f$ Z \f$ a real orthogonal matrix
!> and \f$ E \f$ an antisymmetric real, block-diagonal matrix, with \f$ 2\times 2 \f$ or \f$ 1\times 1 \f$
!> diagonal blocks.
!>
!> Only the lower triangle of \f$ A \f$ is referenced. On output, \p A contains the \f$ Z \f$ matrix,
!> and \p E contains the essential elements of \f$ E \f$, which are pairs of opposite sign
!> for the \f$ 2\times 2 \f$, or zero for the \f$ 1\times 1 \f$ blocks.
!>
!> @param[in]     N    Size of the square matrix
!> @param[in,out] A    Antisymmetric real matrix, it is replaced by its real Schur vectors
!> @param[out]    E    Vector with the real Schur form of \p A
!> @param[out]    ierr Return code (0 if success)
!***********************************************************************

subroutine schur_skew(n,A,E,ierr)
! This code is an adaptation of the code published in:
!   R. C. Ward, L. J. Gray
!   Eigensystem computation for skew-symmetric matrices and a class of symmetric matrices
!   Technical report ORNL/CSD-9 TRN: 76-016851, 1976. Oak Ridge National Laboratory, Tennessee (USA)
!   doi:10.2172/7350249

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(inout) :: A(n,n)
real(kind=wp), intent(out) :: E(n)
integer(kind=iwp), intent(out) :: ierr
integer(kind=iwp) :: i, its, j, k, l, l0, l0m1, lm1, ls, m, m0, mm1
real(kind=wp) :: c, f, g, h, p, q, r, s, sc, test
logical(kind=iwp) :: break, isodd, skip
real(kind=wp), allocatable :: Z(:,:)
integer(kind=iwp), parameter :: maxit = 30
real(kind=wp), parameter :: eps = epsilon(eps)
real(kind=wp), external :: ddot_

ierr = 0

! Reduce to tridiagonal

do i=n,2,-1
  l = i-1
  ! normalize row
  sc = Zero
  do j=1,l
    sc = sc+abs(A(i,j))
  end do
  if (sc < eps) then
    E(i) = Zero
    A(i,i) = Zero
    cycle
  end if
  ! compute elements of U vector
  A(i,1:l) = A(i,1:l)/sc
  ! sum( A(i,1:l) * A(i,1:l) )
  h = ddot_(l,A(i,1),n,A(i,1),n)
  f = A(i,l)
  g = -sign(sqrt(h),f)
  h = h-f*g
  E(i) = sc*g
  A(i,l) = f-g
  A(i,i) = sqrt(h)
  if (l > 1) then
    ! compute elements of A*U/h
    h = One/h
    do j=1,l
      ! sum( A(j,j+1:l) * A(i,j+1:l) ) - sum( A(j+1:l,j) * A(i,j+1:l) )
      E(j) = (ddot_(j-1,A(j,1),n,A(i,1),n)-ddot_(l-j,A(j+1:l,j),1,A(i,j+1),n))*h
    end do
    ! compute reduced A
    do j=2,l
      A(j,1:j-1) = A(j,1:j-1)+A(i,j)*E(1:j-1)-E(j)*A(i,1:j-1)
    end do
  end if
  A(i,1:i) = sc*A(i,1:i)
end do
! enforce conventional value
E(1) = Zero

! Diagonalize skew-symmetric tridiagonal matrix

! place identity matrix in Z
call mma_allocate(Z,n,n,label='Z')
call unitmat(Z,n)

m = n
mm1 = m-1
E(1) = Zero
its = 0

outer: do while (m >= 2)
  m0 = m

  ! search for next submatrix to solve (matrix splitting)
  f = Zero
  break = .false.
  do i=1,mm1
    j = m-i
    g = abs(E(j+1))
    if (g <= eps*(abs(E(j))+f)) then
      break = .true.
      exit
    end if
    f = g
  end do
  if (.not. break) j = 0

  l0 = j+2
  l0m1 = j+1
  if (l0m1 == m) then
    ! iteration converged to one zero eigenvalue
    E(m) = Zero
    m = mm1
    its = 0
    mm1 = m-1
    cycle outer
  end if

  ! place correct sign on identity diagonals
  do i=l0m1,m,4
    Z(i,i) = -Z(i,i)
    if (i+3 > m) exit
    Z(i+3,i+3) = -Z(i+3,i+3)
  end do

  if (l0 == m) then
    ! iteration converged to eigenvalue pair
    E(mm1) = E(m)
    E(m) = -E(m)
    m = m-2
    its = 0
    mm1 = m-1
    cycle outer
  end if

  isodd = mod(m-l0,2) == 1
  l = l0
  if (isodd) then

    ! find zero eigenvalue of odd ordered submatrices
    c = Zero
    s = -One
    do i=l0,mm1,2
      k = mm1+l0-i
      q = -s*E(k+1)
      E(k+1) = c*E(k+1)
      if (abs(E(k)) > abs(q)) then
        s = q/E(k)
        r = sqrt(One+s**2)
        E(k) = E(k)*r
        c = One/r
        s = s*c
      else
        c = E(k)/q
        r = sqrt(c**2+One)
        E(k) = q*r
        s = One/r
        c = c*s
      end if

      ! accumulate transformations for eigenvectors
      Z(k-1,m) = -s*Z(k-1,k-1)
      Z(k-1,k-1) = c*Z(k-1,k-1)
      Z(k+1:m:2,k-1) = s*Z(k+1:m:2,m)
      Z(k+1:m:2,m) = c*Z(k+1:m:2,m)

    end do
    m = mm1
    mm1 = m-1
    if (l0 == m) then
      ! iteration converged to eigenvalue pair
      E(mm1) = E(m)
      E(m) = -E(m)
      m = m-2
      its = 0
      mm1 = m-1
      cycle outer
    end if

  end if

  skip = .not. isodd
  do
    if (.not. skip) then
      ! check for convergence or small subdiagonal element
      break = .false.
      do i=l0,mm1,2
        k = mm1+l0-i
        l = k+1
        test = eps*(abs(E(l))+abs(E(k-1)))
        if (abs(E(k)) <= test) then
          break = .true.
          exit
        end if
      end do
      if (.not. break) l = l0
      if (l == m) then
        ! iteration converged to eigenvalue pair
        do
          E(mm1) = E(m)
          E(m) = -E(m)
          m = m-2
          mm1 = m-1
          if (m /= l0) exit
        end do
        its = 0
        if (m > l0) cycle
        cycle outer
      end if
    end if
    skip = .false.

    ! form shift
    its = its+1
    if (its > maxit) then
      ! error exit
      ierr = m
      exit outer
    end if
    f = E(m-3)
    g = E(m-2)
    c = E(mm1)
    s = E(m)
    p = ((c-f)*(c+f)+(s-g)*(s+g))/(Two*g*c)
    r = sqrt(p*p+One)
    q = (g/(p+sign(r,p)))-c
    f = E(l)
    lm1 = l-1
    E(lm1) = ((f-s)*(f+s)+c*q)/f

    ! perform one implicit QR iteration on Cholesky factor
    ls = l0m1
    c = One
    s = One
    do i=l,mm1
      q = s*E(i+1)
      E(i+1) = c*E(i+1)
      if (abs(E(i-1)) > abs(q)) then
        s = q/E(i-1)
        r = sqrt(One+s**2)
        E(i-1) = E(i-1)*r
        c = One/r
        s = s*c
      else
        c = E(i-1)/q
        r = sqrt(c**2+One)
        E(i-1) = q*r
        s = One/r
        c = c*s
      end if
      f = E(i+1)
      E(i+1) = -s*E(i)+c*f
      E(i) = c*E(i)+s*f

      ! accumulate transformations for eigenvectors
      do j=ls,m0,2
        f = Z(j,i+1)
        Z(j,i+1) = -s*Z(j,i-1)+c*f
        Z(j,i-1) = c*Z(j,i-1)+s*f
      end do
      if (ls == l0m1) then
        ls = l0
      else
        ls = l0m1
      end if
    end do
    E(lm1) = Zero
  end do

end do outer

! Back-transform to full matrix

do i=2,n
  if (A(i,i) == Zero) cycle
  do j=1,n
    ! sum( A(i,1:i-1) * Z(1:i-1,j) ) / A(i,i)**2
    s = (ddot_(i-1,A(i,1),n,Z(1:i-1,j),1)/A(i,i))/A(i,i)
    Z(1:i-1,j) = Z(1:i-1,j)-s*A(i,1:i-1)
  end do
end do
A(:,:) = Z

call mma_deallocate(Z)

return

end subroutine schur_skew
