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

subroutine MINV_INNER(ARRAY,ARRINV,ISING,DET,NDIM,A,BUF,B,IPIV,JPIV)
!subroutine Dool(NDIM,MDIM,N,M,A,B,DET,IPIV,JPIV,BUF)
! SOLVES THE MATRIX EQUATION AX=B BY DOOLITTLE'S METHOD
! ACTUAL DIMENSIONS ARE N*N AND N*M
! ALLOCATED DIMENSIONS ARE NDIM*NDIM AND NDIM*MDIM
! A AND B ARE DESTROYED, AND X IS RETURNED AS MATRIX B
!                                   (MALMQUIST 82-11-12)
!                                    (UPDATE 83-09-28)

implicit real*8(a-h,o-z)
#include "real.fh"
! --- global variables ---
real*8 ARRAY(NDIM,NDIM), ARRINV(NDIM,NDIM)
real*8 A(NDIM,NDIM), BUF(NDIM), B(NDIM,NDIM)
integer IPIV(NDIM), JPIV(NDIM)

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

do I=1,NDIM
  do J=1,NDIM
    A(I,J) = ARRAY(I,J)
    B(I,J) = Zero
  end do
  B(I,I) = One
end do

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
      if (AMAX > AM) Go To 20
      AMAX = AM
      IP = K
      JP = L
20    continue
    end do
  end do
  if (IP == I) Go To 3
  DET = -DET
  IDUM = IPIV(I)
  IPIV(I) = IPIV(IP)
  IPIV(IP) = IDUM
3 if (JP == I) Go To 4
  DET = -DET
  IDUM = JPIV(I)
  JPIV(I) = JPIV(JP)
  JPIV(JP) = IDUM
4 IP = IPIV(I)
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
    SUM = B(IP,J)
    do K=1,I-1
      SUM = SUM-A(IP,JPIV(K))*B(IPIV(K),J)
    end do
    B(IP,J) = SUM
  end do
end do

! SECOND RESUBSTITUTION STEP:

do J=1,M
  do I=N,1,-1
    IP = IPIV(I)
    SUM = B(IP,J)
    do K=I+1,N
      SUM = SUM-A(IP,JPIV(K))*B(IPIV(K),J)
    end do
    if (BUF(I) /= Zero) SUM = SUM/BUF(I)
    B(IP,J) = SUM
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

do I=1,NDIM
  do J=1,NDIM
    ARRINV(I,J) = B(I,J)
  end do
end do

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(ISING)

end subroutine MINV_INNER
