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

subroutine MINV(ARRAY,ARRINV,DET,NDIM)
!subroutine Dool(NDIM,MDIM,N,M,A,B,DET,IPIV,JPIV,BUF)
! SOLVES THE MATRIX EQUATION AX=B BY DOOLITTLE'S METHOD
! ACTUAL DIMENSIONS ARE N*N AND N*M
! ALLOCATED DIMENSIONS ARE NDIM*NDIM AND NDIM*MDIM
! A AND B ARE DESTROYED, AND X IS RETURNED AS MATRIX B
!                                   (MALMQUIST 82-11-12)
!                                    (UPDATE 83-09-28)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NDIM
real(kind=wp), intent(in) :: ARRAY(NDIM,NDIM)
real(kind=wp), intent(out) :: ARRINV(NDIM,NDIM), DET
integer(kind=iwp) :: I, IDUM, ip, J, jp, K, KP, L, LP, M, N
real(kind=wp) :: AM, AMAX, C, DIAG, RSUM
integer(kind=iwp), allocatable :: IPIV(:), JPIV(:)
real(kind=wp), allocatable :: A(:,:), B(:,:), BUF(:)

call mma_allocate(A,NDIM,NDIM,label='A')
call mma_allocate(B,NDIM,NDIM,label='B')
call mma_allocate(BUF,NDIM,label='BUF')
call mma_allocate(IPIV,NDIM,label='IPIV')
call mma_allocate(JPIV,NDIM,label='JPIV')

! EQUATION IS SOLVED BY FACTORIZING A=L*R IN SAME SPACE AS A.
! PIVOTING IS ACHIEVED BY INDIRECT INDEXING.
! FIRST PIVOTING INDICES ARE ASSIGNED START VALUES.

! set N and M to NDIM since this subroutine is modified to only
! deal with square matrices.

N = NDIM
M = NDIM
!MDIM = NDIM

! Move ARRAY to allocation A and set the B matrix equal to the
! unity matrix

A(:,:) = ARRAY
call unitmat(B,NDIM)

! Let's go!!

ip = -1
jp = -1
do I=1,N
  IPIV(I) = I
  JPIV(I) = I
end do
DET = One
do I=1,N

  ! NOW FIND BETTER PIVOT ELEMENT:

  AMAX = -One
  do K=I,N
    do L=I,N
      AM = abs(A(IPIV(K),JPIV(L)))
      if (AMAX > AM) cycle
      AMAX = AM
      IP = K
      JP = L
    end do
  end do
  if (IP /= I) then
    DET = -DET
    IDUM = IPIV(I)
    IPIV(I) = IPIV(IP)
    IPIV(IP) = IDUM
  end if
  if (JP /= I) then
    DET = -DET
    IDUM = JPIV(I)
    JPIV(I) = JPIV(JP)
    JPIV(JP) = IDUM
  end if
  IP = IPIV(I)
  JP = JPIV(I)
  DIAG = A(IP,JP)
  BUF(I) = DIAG
  DET = DET*DIAG
  do K=I+1,N
    KP = IPIV(K)
    C = A(KP,JP)
    if (DIAG /= Zero) C = C/DIAG
    A(KP,JP) = C
    do L=I+1,N
      LP = JPIV(L)
      A(KP,LP) = A(KP,LP)-C*A(IP,LP)
    end do
  end do
end do

! FIRST RESUBSTITUTION STEP:

do J=1,M
  do I=2,N
    IP = IPIV(I)
    RSUM = B(IP,J)
    do K=1,I-1
      RSUM = RSUM-A(IP,JPIV(K))*B(IPIV(K),J)
    end do
    B(IP,J) = RSUM
  end do
end do

! SECOND RESUBSTITUTION STEP:

do J=1,M
  do I=N,1,-1
    IP = IPIV(I)
    RSUM = B(IP,J)
    do K=I+1,N
      RSUM = RSUM-A(IP,JPIV(K))*B(IPIV(K),J)
    end do
    if (BUF(I) /= Zero) RSUM = RSUM/BUF(I)
    B(IP,J) = RSUM
  end do
end do

! REORGANIZATION PART:

do J=1,M
  do I=1,N
    BUF(I) = B(IPIV(I),J)
  end do
  do I=1,N
    B(JPIV(I),J) = BUF(I)
  end do
end do

! Move the result to location ARRINV

ARRINV(:,:) = B

call mma_deallocate(A)
call mma_deallocate(B)
call mma_deallocate(BUF)
call mma_deallocate(IPIV)
call mma_deallocate(JPIV)

return

end subroutine MINV
