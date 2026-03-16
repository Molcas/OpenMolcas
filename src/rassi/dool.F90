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
! Copyright (C) 1982,1983, Per Ake Malmqvist                           *
!***********************************************************************

! SOLVES THE MATRIX EQUATION AX=B BY DOOLITTLE'S METHOD
! ACTUAL DIMENSIONS ARE N*N AND N*M
! ALLOCATED DIMENSIONS ARE NDIMEN*NDIMEN AND NDIMEN*MDIM
! A AND B ARE DESTROYED, AND X IS RETURNED AS MATRIX B
!                                   (MALMQUIST 82-11-12)
!                                    (UPDATE 83-09-28)
subroutine DOOL(NDIMEN,MDIM,N,M,A,B,DET,IPIV,JPIV,BUF)

use definitions, only: iwp, wp
use constants, only: One

implicit none
integer(kind=iwp), intent(in) :: NDIMEN, MDIM, N, M
real(kind=wp), intent(inout) :: A(NDIMEN,NDIMEN), B(NDIMEN,MDIM)
real(kind=wp), intent(out) :: DET
integer(kind=iwp), intent(out) :: IPIV(NDIMEN), JPIV(NDIMEN)
real(kind=wp), intent(out) :: BUF(NDIMEN)
real(kind=wp) AM, AMAX, C, DIAG, SUMMA
integer(kind=iwp) I, IDUM, IP, J, JP, K, KP, L, LP

! EQUATION IS SOLVED BY FACTORIZING A=L*R IN SAME SPACE AS A.
! PIVOTING IS ACHIEVED BY INDIRECT INDEXING.
! FIRST PIVOTING INDICES ARE ASSIGNED START VALUES.

IP = 0 ! dummy initialize
JP = 0 ! dummy initialize
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
      if (AMAX <= AM) then
        AMAX = AM
        IP = K
        JP = L
      end if
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
    C = A(KP,JP)/DIAG
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
    SUMMA = B(IP,J)
    do K=1,I-1
      SUMMA = SUMMA-A(IP,JPIV(K))*B(IPIV(K),J)
    end do
    B(IP,J) = SUMMA
  end do
end do

! SECOND RESUBSTITUTION STEP:

do J=1,M
  do I=N,1,-1
    IP = IPIV(I)
    SUMMA = B(IP,J)
    do K=I+1,N
      SUMMA = SUMMA-A(IP,JPIV(K))*B(IPIV(K),J)
    end do
    B(IP,J) = SUMMA/BUF(I)
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

end subroutine DOOL
