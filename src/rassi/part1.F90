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
! Copyright (C) 1984,1989, Per Ake Malmqvist                           *
!***********************************************************************

subroutine PART1(NDIMEN,NBLOCK,NSIZE,SXY,B,A,SCR,IPIV,BUF)
! PURPOSE: SEE SUBROUTINE PART.
! SUBDIVISION INTO TWO LEVELS OF ROUTINE CALLS IS MERELY TO
! FACILITATE HANDLING OF SYMMETRY AND INDEXING.
! ORIGINAL VERSION, MALMQUIST 84-04-04
! RASSCF VERSION,   MALMQUIST 89-11-15

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NDIMEN, NBLOCK, NSIZE(NBLOCK)
real(kind=wp), intent(in) :: SXY(NDIMEN,NDIMEN)
real(kind=wp), intent(out) :: B(NDIMEN,NDIMEN), A(NDIMEN,NDIMEN), BUF(NDIMEN)
real(kind=wp), intent(inout) :: SCR(NDIMEN,NDIMEN)
integer(kind=iwp), intent(out) :: IPIV(NDIMEN,2)
integer(kind=iwp) :: I, J, K, L, LIM1, LIM2, LIM3, NSZ
real(kind=wp) :: DET, T

! INITIALIZE A = INVERSE OF SXY, AND B = UNIT MATRIX:
SCR(:,:) = SXY(:,:)
call unitmat(A,NDIMEN)
call unitmat(B,NDIMEN)
call DOOL(NDIMEN,NDIMEN,NDIMEN,NDIMEN,SCR,A,DET,IPIV(:,1),IPIV(:,2),BUF)
!---------------------------------------------------------------
! LOOP BACKWARDS OVER THE BLOCKS. KEEP THREE LIMITS UPDATED:
! LIM1= POS IMMEDIATELY BEFORE CURRENT BLOCK, LIM2= BEGINNING
! OF CURRENT BLOCK, AND LIM3=END OF CURRENT BLOCK.
LIM1 = NDIMEN
do K=NBLOCK,2,-1
  NSZ = NSIZE(K)
  LIM3 = LIM1
  LIM1 = LIM1-NSZ
  LIM2 = LIM1+1
  !---------------------------------------------------------------
  ! CALCULATE (INVERSE OF CURRENT A-BLOCK)*(ALL TO THE LEFT OF IT)
  ! AND PUT IT INTO B-MATRIX. THEN CLEAR ALL TO THE LEFT OF THE A-BLOCK.
  SCR(LIM2:LIM3,LIM2:LIM3) = A(LIM2:LIM3,LIM2:LIM3)
  B(LIM2:LIM3,1:LIM1) = A(LIM2:LIM3,1:LIM1)
  A(LIM2:LIM3,1:LIM1) = Zero
  call DOOL(NDIMEN,NDIMEN,NSZ,LIM1,SCR(LIM2,LIM2),B(LIM2,1),DET,IPIV(1,1),IPIV(1,2),BUF)
  !---------------------------------------------------------------
  ! THEN UPDATE THE COLUMNS OF A TO THE LEFT OF THE CURRENT BLOCK:
  do J=1,LIM1
    do I=1,LIM1
      A(I,J) = A(I,J)-sum(B(LIM2:LIM3,J)*A(I,LIM2:LIM3))
    end do
  end do
end do
! TRANSPOSE MATRIX B:
do I=1,NDIMEN-1
  do J=I,NDIMEN
    T = B(I,J)
    B(I,J) = B(J,I)
    B(J,I) = T
  end do
end do
!---------------------------------------------------------------
! COMBINED LU-PARTITIONING AND UNITARY TRANSFORMATION OF A AND B:
call LU2(NDIMEN,NBLOCK,NSIZE,B,A,BUF)

! LU PARTITIONING OF THE MATRICES WAS DONE IN-PLACE.
! NOW CHANGE SIGN OF LOWER-TRIANGULAR PARTS AND
! INVERT UPPER-TRIANGULAR PARTS, AS INDICATED IN (VI.6):

do J=1,NDIMEN-1
  A(J+1:NDIMEN,J) = -A(J+1:NDIMEN,J)
  B(J+1:NDIMEN,J) = -B(J+1:NDIMEN,J)
end do
do L=NDIMEN,1,-1
  A(L,L) = One/A(L,L)
  B(L,L) = One/B(L,L)
  A(L,L+1:NDIMEN) = A(L,L)*A(L,L+1:NDIMEN)
  B(L,L+1:NDIMEN) = B(L,L)*B(L,L+1:NDIMEN)
  do K=1,L-1
    A(K,L+1:NDIMEN) = A(K,L+1:NDIMEN)-A(K,L)*A(L,L+1:NDIMEN)
    B(K,L+1:NDIMEN) = B(K,L+1:NDIMEN)-B(K,L)*B(L,L+1:NDIMEN)
    A(K,L) = -A(L,L)*A(K,L)
    B(K,L) = -B(L,L)*B(K,L)
  end do
end do

end subroutine PART1
