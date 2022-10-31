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

subroutine Givens(H,U,n,nv)
!***********************************************************************
!                                                                      *
! This routine transforms a symmetric matrix to tridiagonal form using *
! Givens rotations. The matrix is stored in lower triangular form.     *
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
!----------------------------------------------------------------------*
! Parameters                                                           *
! zThr - Threshold for when an element is regarded as zero             *
!----------------------------------------------------------------------*

use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: n, nv
real(kind=wp) :: H(*), U(nv,n)
integer(kind=iwp) :: i, ii, ij, ik, im, iSkip, j, jj, jk, jm, k, m
real(kind=wp) :: Hii, Hij, Hjj, p, q, tmp
real(kind=wp), parameter :: zThr = 1.0e-16_wp

!----------------------------------------------------------------------*
!                                                                      *
! i,j - rotate around these indices                                    *
! i,k - eliminate this element (k=j-1)                                 *
!                                                                      *
!----------------------------------------------------------------------*
do j=2,n-1
  do i=j+1,n
    k = j-1
    ii = i*(i+1)/2
    jj = j*(j+1)/2
    ij = j+i*(i-1)/2
    ik = k+i*(i-1)/2
    jk = k+j*(j-1)/2
    Hii = H(ii)
    Hjj = H(jj)
    Hij = H(ij)

    iSkip = 0
    if (abs(H(ik)) < zThr) then
      p = One
      q = Zero
      iSkip = 1
    else if (abs(H(jk)) < zThr) then
      p = Zero
      q = One
    else if (abs(H(jk)) < abs(H(ik))) then
      tmp = H(jk)/H(ik)
      p = tmp/sqrt(One+tmp*tmp)
      q = sqrt(One-p*p)
      if (p < Zero) then
        p = -p
        q = -q
      end if
    else
      tmp = H(ik)/H(jk)
      q = tmp/sqrt(One+tmp*tmp)
      p = sqrt(One-q*q)
    end if
    if (iSkip == 1) cycle

    do m=1,n
      if (m < j) then
        im = m+i*(i-1)/2
        jm = m+j*(j-1)/2
      else if (m < i) then
        im = m+i*(i-1)/2
        jm = j+m*(m-1)/2
      else
        im = i+m*(m-1)/2
        jm = j+m*(m-1)/2
      end if
      tmp = p*H(im)-q*H(jm)
      H(jm) = q*H(im)+p*H(jm)
      H(im) = tmp
    end do

    H(ii) = p*p*Hii+q*q*Hjj-Two*p*q*Hij
    H(jj) = q*q*Hii+p*p*Hjj+Two*p*q*Hij
    H(ij) = (p*p-q*q)*Hij+p*q*(Hii-Hjj)
    H(ik) = Zero
    do m=1,nv
      tmp = p*U(m,i)-q*U(m,j)
      U(m,j) = q*U(m,i)+p*U(m,j)
      U(m,i) = tmp
    end do
  end do
end do

!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine Givens
