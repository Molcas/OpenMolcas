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
! Copyright (C) 1986,1995, Bernd Artur Hess                            *
!               2005, Jesper Wisborg Krogh                             *
!***********************************************************************

!#NUMPAC#MINVD               REVISED ON 1984-11-30
subroutine MINVD(A,KA,N,EPS,ILL)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: KA, N
real(kind=wp), intent(inout) :: A(KA,N)
real(kind=wp), intent(in) :: EPS
integer(kind=iwp), intent(out) :: ILL
integer(kind=iwp) :: I, IM1, IM2, J, JM1, JP1, K, M, NM1
real(kind=wp) :: AA, AM, P, S, W
integer(kind=iwp), allocatable :: MX(:)

if ((N < 1) .or. (N > KA) .or. (EPS <= Zero)) then
  ILL = 30000
  return
end if

call mma_allocate(MX,N,label='MX')

! LU DECOMPOSITION
JM1 = 0 ! dummy initialize
M = 0 ! dummy initialize
NM1 = N-1
do J=1,N
  if (J /= 1) then
    JM1 = J-1
    do I=1,JM1
      M = MX(I)
      S = A(M,J)
      A(M,J) = A(I,J)
      if (I /= 1) then
        IM1 = I-1
        do K=1,IM1
          S = A(I,K)*A(K,J)+S
        end do
      end if
      A(I,J) = S
    end do
  end if
  AM = Zero
  do I=J,N
    S = A(I,J)
    if (J /= 1) then
      do K=1,JM1
        S = A(I,K)*A(K,J)+S
      end do
      A(I,J) = S
    end if
    AA = abs(S)
    if (AA > AM) then
      AM = AA
      M = I
    end if
  end do
  if (AM < EPS) then
    ILL = J
    call mma_deallocate(MX)
    return
  end if
  MX(J) = M
  if (M /= J) then
    do K=1,J
      W = A(M,K)
      A(M,K) = A(J,K)
      A(J,K) = W
    end do
  end if
  if (J /= N) then
    JP1 = J+1
    W = -A(J,J)
    do I=JP1,N
      A(I,J) = A(I,J)/W
    end do
  end if
end do
if (N > 2) then
  ! INPLACE INVERSION OF L-COMPONENT
  do I=3,N
    IM1 = I-1
    IM2 = I-2
    do J=1,IM2
      S = A(I,J)
      JP1 = J+1
      do K=JP1,IM1
        S = A(I,K)*A(K,J)+S
      end do
      A(I,J) = S
    end do
  end do
end if
! INPLACE INVERSION OF U-COMPONENT
A(1,1) = One/A(1,1)
if (N == 1) then
  ILL = 0
  call mma_deallocate(MX)
  return
end if
do J=2,N
  A(J,J) = One/A(J,J)
  P = -A(J,J)
  JM1 = J-1
  do I=1,JM1
    S = Zero
    do K=I,JM1
      S = A(I,K)*A(K,J)+S
    end do
    A(I,J) = S*P
  end do
end do
! INPLACE MULTIPLICATION OF L AND U COMPONENT
do J=1,NM1
  JP1 = J+1
  do I=1,J
    S = A(I,J)
    do K=JP1,N
      S = A(I,K)*A(K,J)+S
    end do
    A(I,J) = S
  end do
  do I=JP1,N
    S = Zero
    do K=I,N
      S = A(I,K)*A(K,J)+S
    end do
    A(I,J) = S
  end do
end do
! INTERCHANGE OF COLUMNS
J = NM1
do
  M = MX(J)
  if (M /= J) then
    do I=1,N
      W = A(I,M)
      A(I,M) = A(I,J)
      A(I,J) = W
    end do
  end if
  J = J-1
  if (J < 1) exit
end do

call mma_deallocate(MX)

return

end subroutine MINVD
