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
! Copyright (C) 1993, Markus P. Fuelscher                              *
!***********************************************************************

subroutine Schmidt(N,S,C,Temp,M)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Schmidt orthogonalization such that the latest vector (in chro-  *
!     nological order) remains unchanged and the previous are ortho-   *
!     gonalized relative to it.                                        *
!                                                                      *
!     calling arguments:                                               *
!     N       : Type integer, input.                                   *
!               Dimensions of the overlap matrix.                      *
!     S       : Type double precision real, input.                     *
!               Overlap matrix.                                        *
!     C       : Type double precision real, output.                    *
!               Matrix containing the transformation vectors.          *
!     Temp    : Type double precision real, input/output.              *
!               Scratch area.                                          *
!     M       : Type integer, output.                                  *
!               This is the number of linearly independent basis       *
!               vectors that span space the metric given by the        *
!               overlap matrix.                                        *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher, University of Lund, Sweden, 1993                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N
integer(kind=iwp), intent(out) :: M
real(kind=wp), intent(in) :: S(N,N)
real(kind=wp), intent(out) :: C(N,N), Temp(N)
integer(kind=iwp) :: i, j, k
real(kind=wp) :: Alpha, rSum
logical(kind=iwp), parameter :: forward = .true.

M = 0
C(:,:) = Zero
do i=1,N
  C(i,i) = One/sqrt(S(i,i))
end do

! forward orthonormalization
! (vectors 1 remains unchanged)
if (forward) then
  do i=1,N
    Alpha = C(i,i)
    do j=1,N
      Temp(j) = S(j,i)*Alpha
    end do
    do j=1,i-1
      rSum = Zero
      do k=1,i
        rSum = rSum+C(k,j)*Temp(k)
      end do
      do k=1,i
        C(k,i) = C(k,i)-rSum*C(k,j)
      end do
    end do
    rSum = Zero
    do k=1,i
      rSum = rSum+C(k,i)*Temp(k)
    end do
    if (rSum > 1.0e-9_wp) then
      M = M+1
      Alpha = One/sqrt(rSum)
      do k=1,i
        C(k,i) = C(k,i)*Alpha
      end do
    else
      do k=1,i
        C(k,i) = Zero
      end do
    end if
  end do
end if

! backward orthonormalization
! (vectors N remains unchanged)
if (.not. forward) then
  do i=N,1,-1
    Alpha = C(i,i)
    do j=1,N
      Temp(j) = S(j,i)*Alpha
    end do
    do j=N,i+1,-1
      rSum = Zero
      do k=j,N
        rSum = rSum+C(k,j)*Temp(k)
      end do
      do k=j,N
        C(k,i) = C(k,i)-rSum*C(k,j)
      end do
    end do
    rSum = Zero
    do k=i,N
      rSum = rSum+C(k,i)*Temp(k)
    end do
    if (rSum > 1.0e-9_wp) then
      M = M+1
      Alpha = One/sqrt(rSum)
      do k=i,N
        C(k,i) = C(k,i)*Alpha
      end do
    else
      do k=i,N
        C(k,i) = Zero
      end do
    end if
  end do
end if

return

end subroutine Schmidt
