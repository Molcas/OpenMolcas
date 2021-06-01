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
! Copyright (C) 2000, Markus P. Fuelscher                              *
!***********************************************************************

subroutine Eigen_Molcas(N,X,D,E)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Compute eigenvalues and eigenvectors of a symmetric real matrix  *
!     by the method of Householder (QR algorithm)                      *
!                                                                      *
!     reference:                                                       *
!     R.S. Martin, C. Reinsch and J.H. Wilkinson                       *
!     Num. Mat. Vol 11, p 181.-195 (1968)                              *
!                                                                      *
!     calling arguments:                                               *
!     N       : Type integer, input.                                   *
!               Dimensions of the matrix X and vectors D and E         *
!     X       : Type real*8 real, input/output                         *
!               on input it is the matrix to diagonalized              *
!               on output it contains the eigenvectors                 *
!     D       : Type real*8 real, output.                              *
!               vector of eigenvalues                                  *
!     E       : Type real*8 real, input/output.                        *
!               Scratch area of length N                               *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher, University of Lund, Sweden, 2000                 *
!     (the subrotine is based on a old implementation written          *
!      by the comp. center in Munich)                                  *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(inout) :: X(N,N), E(N)
real(kind=wp), intent(out) :: D(N)
integer(kind=iwp) :: i, ii, j, k, l
real(kind=wp) :: B, C, dsum, F, G, H, P, R, S, scal, tmp
real(kind=wp), parameter :: eps = 3.0e-17_wp, tol = 1.0e-22_wp

if (N == 1) then
  D(1) = X(1,1)
  X(1,1) = One
  return
end if

! HOUSEHOLDER REDUCTION

do ii=2,N
  i = N+2-ii
  l = i-2
  H = Zero
  G = X(i,i-1)
  if (l > 0) then
    do k=1,l
      H = H+X(i,k)*X(i,k)
    end do
    S = H+G*G
    if (S < Tol) H = Zero
    if (H > Zero) then
      l = l+1
      F = G
      G = sqrt(S)
      if (G > Zero) G = -G
      H = S-F*G
      X(i,i-1) = F-G
      F = Zero
      do j=1,l
        X(j,i) = X(i,j)/H
        dsum = Zero
        do k=1,j
          dsum = dsum+X(j,k)*X(i,k)
        end do
        do k=j+1,l
          dsum = dsum+X(k,j)*X(i,k)
        end do
        E(j) = dsum/H
        F = F+dsum*X(j,i)
      end do
      F = F/(H+H)
      do j=1,l
        E(j) = E(j)-F*X(i,j)
      end do
      do j=1,l
        F = X(i,j)
        scal = E(j)
        do k=1,j
          X(j,k) = X(j,k)-F*E(k)-X(i,k)*scal
        end do
      end do
    end if
  end if
  D(i) = H
  E(i-1) = G
end do

! ACCUMULATION OF TRANSFORMATION MATRICES

D(1) = X(1,1)
X(1,1) = One
do i=2,N
  l = i-1
  if (D(i) > Zero) then
    do j=1,l
      S = Zero
      do k=1,l
        S = S+X(i,k)*X(k,j)
      end do
      do k=1,l
        X(k,j) = X(k,j)-S*X(k,i)
      end do
    end do
  end if
  D(i) = X(i,i)
  X(i,i) = One
  do j=1,l
    X(i,j) = Zero
    X(j,i) = Zero
  end do
end do

! DIAGONALIZATION OF THE TRIDIAGONAL MATRIX

B = Zero
F = Zero
E(N) = Zero
outer: do l=1,N
  H = eps*(abs(D(l))+abs(E(l)))
  if (H > B) B = H
  do j=l,N
    if (abs(E(j)) <= B) then
      if (j == l) then
        D(l) = D(l)+F
        cycle outer
      end if
    end if
  end do
  j = N
  do while (abs(E(l)) > B)
    P = (D(l+1)-D(l))*Half/E(l)
    R = sqrt(P*P+One)
    if (P >= Zero) then
      P = P+R
    else
      P = P-R
    end if
    H = D(l)-E(l)/P
    do i=l,N
      D(i) = D(i)-H
    end do
    F = F+H
    P = D(j)
    C = One
    S = Zero
    do ii=l,j-1
      i = l+j-1-ii
      G = C*E(i)
      H = C*P
      if (abs(P) < abs(E(i))) then
        C = P/E(i)
        R = sqrt(C*C+One)
        E(i+1) = S*E(i)*R
        S = One/R
        C = C/R
      else
        C = E(i)/P
        R = sqrt(C*C+One)
        E(i+1) = S*P*R
        S = C/R
        C = One/R
      end if
      P = C*D(i)-S*G
      D(i+1) = H+S*(C*G+S*D(i))
      do k=1,N
        H = X(k,i+1)
        X(k,i+1) = X(k,i)*S+H*C
        X(k,i) = X(k,i)*C-H*S
      end do
    end do
    E(l) = S*P
    D(l) = C*P
  end do
  D(l) = D(l)+F
end do outer

! ORDERING OF EIGENVALUES

do i=1,N-1
  do j=i+1,N
    if (D(j) < D(i)) then
      tmp = D(j)
      D(j) = D(i)
      D(i) = tmp
      do k=1,N
        tmp = X(k,j)
        X(k,j) = X(k,i)
        X(k,i) = tmp
      end do
    end if
  end do
end do

! FIXING OF SIGN

do i=1,N
  k = 1
  tmp = abs(X(k,i))
  do j=2,N
    if (tmp <= abs(X(j,i))) then
      tmp = abs(X(j,i))
      k = j
    end if
  end do
  if (X(k,i) < Zero) then
    do j=1,N
      X(j,i) = -X(j,i)
    end do
  end if
end do

return

end subroutine Eigen_Molcas
