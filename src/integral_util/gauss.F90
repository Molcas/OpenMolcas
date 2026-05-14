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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!***********************************************************************

subroutine Gauss(n,lDim,A,X,C)
!***********************************************************************
!                                                                      *
!     purpose: Solve a set of linear equations using Gauss method      *
!                            A*X = C                                   *
!                                                                      *
!     input:                                                           *
!       A,C     : input matrices                                       *
!       n       : size of the set of equations                         *
!       lDim    : leading dimension of A                               *
!                                                                      *
!     output:                                                          *
!       X       : solutions                                            *
!                                                                      *
!     called from: Diis, MinDns                                        *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
!     University of Lund, Sweden, 1992                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n, lDim
real(kind=wp), intent(inout) :: A(lDim,n)
real(kind=wp), intent(out) :: X(n)
real(kind=wp), intent(in) :: C(n)
integer(kind=iwp) :: i, j, k
real(kind=wp) :: Fact
real(kind=wp), allocatable :: Swap(:)

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
call mma_allocate(Swap,n,label='Swap')

X(:) = C(:)
do i=1,n-1
  k = i
  do j=i+1,n
    if (abs(A(k,i)) < abs(A(j,i))) k = j
  end do
  if (k /= i) then
    !write(u6,'(A,2I3)') ' Swapping:',i,k
    Swap(i:) = A(i,i:)
    A(i,i:) = A(k,i:)
    A(k,i:) = Swap(i:)
    Swap(1) = X(i)
    X(i) = X(k)
    X(k) = Swap(1)
  end if
  do k=i+1,n
    Fact = A(k,i)/A(i,i)
    A(k,i+1:) = A(k,i+1:)-Fact*A(i,i+1:)
    X(k) = X(k)-Fact*X(i)
  end do
end do
X(n) = X(n)/A(n,n)
do i=n-1,1,-1
  X(i) = (X(i)-sum(A(i,i+1:)*X(i+1:)))/A(i,i)
end do

call mma_deallocate(Swap)

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

return

end subroutine Gauss
