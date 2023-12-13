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

subroutine MATINV(A,B,N,L,I_DIM)
! IF L=0 RETURNS INVERSE OF A IN A, IF L=1 SOLUTION OF AX=B IN B

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: N, L, I_DIM
real(kind=wp), intent(inout) :: A(I_DIM,I_DIM), B(I_DIM)
integer(kind=iwp) :: I, IC, IR, J, K
real(kind=wp) :: AMAX, D
integer(kind=iwp), allocatable :: I_N(:,:), IP(:)

call mma_allocate(I_N,I_DIM,2,label='I_N')
call mma_allocate(IP,I_DIM,label='IP')

IR = 0
IC = 0
D = One
do I=1,N
  IP(I) = 0
end do
do I=1,N
  AMAX = Zero
  do J=1,N
    if (IP(J) > 0) cycle
    if (IP(J) < 0) then
      write(u6,105)
      call Abend()
    end if
    do K=1,N
      if (IP(K) == 1) cycle
      if (IP(K) > 1) then
        write(u6,105)
        call Abend()
      end if
      if (abs(A(J,K)) <= AMAX) cycle
      IR = J
      IC = K
      AMAX = abs(A(J,K))
    end do
  end do
  IP(IC) = IP(IC)+1
  if (AMAX <= 1.0e-30_wp) then
    write(u6,105)
    call Abend()
  end if
  if (IR /= IC) then
    D = -D
    do K=1,N
      AMAX = A(IR,K)
      A(IR,K) = A(IC,K)
      A(IC,K) = AMAX
    end do
    if (L /= 0) then
      AMAX = B(IR)
      B(IR) = B(IC)
      B(IC) = AMAX
    end if
  end if
  I_N(I,1) = IR
  I_N(I,2) = IC
  AMAX = A(IC,IC)
  D = D*AMAX
  A(IC,IC) = One
  do K=1,N
    A(IC,K) = A(IC,K)/AMAX
  end do
  if (L /= 0) B(IC) = B(IC)/AMAX
  do J=1,N
    if (J == IC) cycle
    AMAX = A(J,IC)
    A(J,IC) = Zero
    do K=1,N
      A(J,K) = A(J,K)-A(IC,K)*AMAX
    end do
    if (L /= 0) B(J) = B(J)-B(IC)*AMAX
  end do
end do
if (L /= 1) then
  do I=1,N
    J = N+1-I
    if (I_N(J,1) == I_N(J,2)) cycle
    IR = I_N(J,1)
    IC = I_N(J,2)
    do K=1,N
      AMAX = A(K,IR)
      A(K,IR) = A(K,IC)
      A(K,IC) = AMAX
    end do
  end do
end if

call mma_deallocate(I_N)
call mma_deallocate(IP)

return

105 format(' * '/' *  SINGULAR MATRIX')

end subroutine MATINV
