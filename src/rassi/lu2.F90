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
! Copyright (C) 1984, Per Ake Malmqvist                                *
!***********************************************************************

subroutine LU2(NDIMEN,NBLOCK,NSIZE,CXA,CYB,SCR)
! GIVES A SIMULTANEOUS LU-PARTITIONING OF MATRICES CXA,CYB IN THE
! SENSE THAT CXA*X = L1*U1 AND CYB*X = L2*U2, WHERE X IS A BLOCK
! UNITARY MATRIX. X IS NEVER FORMED, BUT AT EACH STEP SUCH A
! TRANSFORMATION IS APPLIED THAT THE LU PARTITIONING CAN BE
! APPLIED IN AN OPTIMAL WAY -- THIS IS REFERRED TO AS A
! UNITARY PSEUDO-PIVOTATION. MATRICES CXA AND CYB ARE DESTROYED,
! AND WILL CONTAIN THE NONTRIVIAL ELEMENTS OF THE TRIANGULAR
! MATRICES.
!                                        ( MALMQUIST 84-01-16 )

use Constants, only: One, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NDIMEN, NBLOCK, NSIZE(NBLOCK)
real(kind=wp), intent(inout) :: CXA(NDIMEN,NDIMEN), CYB(NDIMEN,NDIMEN), SCR(NDIMEN)
integer(kind=iwp) :: I, IBLOCK, IEND, II, ISTA, NTOT
real(kind=wp) :: S, S1, S2, S3, X, X1, X2
real(kind=wp), parameter :: THR = 1.0e-6_wp

NTOT = sum(NSIZE(:))
if (NTOT /= NDIMEN) write(u6,*) ' ERROR: NTOT /= NDIMEN IN LU2.'
IEND = 0
do IBLOCK=1,NBLOCK
  ISTA = IEND+1
  IEND = IEND+NSIZE(IBLOCK)
  do II=ISTA,IEND
    S1 = sum(CXA(II,II:IEND)**2)
    S2 = sum(CYB(II,II:IEND)**2)
    S3 = sum(CXA(II,II:IEND)*CYB(II,II:IEND))
    if ((S1 < THR) .or. (S2 < THR)) then
      write(u6,*) ' RASSI CANNOT CONTINUE. THE PROBLEM AT HAND'
      write(u6,*) ' IS PROBABLY NOT SOLUBLE. THE TWO ORBITAL'
      write(u6,*) ' SPACES ARE TOO DISSIMILAR.'
      write(u6,*) ' LU PARTITIONING IS TROUBLESOME. DIAGONAL'
      write(u6,*) ' ELEMENT NR.',II,' IS TOO SMALL:'
      write(u6,*) ' IN MATRIX CXA IT IS',sqrt(S1)
      write(u6,*) ' IN MATRIX CYB IT IS',sqrt(S2)
      write(u6,*) ' EVEN AFTER OPTIMAL PIVOT-TRANSFORMATION.'
      call ABEND()
    end if
    X1 = One/sqrt(S1)
    X2 = One/sign(sqrt(S2),S3)
    SCR(II:IEND) = X1*CXA(II,II:IEND)+X2*CYB(II,II:IEND)
    S = Two*(One+X1*X2*S3)
    X = One/sign(sqrt(S),SCR(II))
    SCR(II:IEND) = X*SCR(II:IEND)
    X = One/(One+SCR(II))
    do I=1,IEND
      S2 = sum(CXA(I,II+1:IEND)*SCR(II+1:IEND))
      S = S2+CXA(I,II)*SCR(II)
      S2 = CXA(I,II)+X*S2
      CXA(I,II) = S
      CXA(I,II+1:IEND) = CXA(I,II+1:IEND)-S2*SCR(II+1:IEND)
    end do
    do I=1,IEND
      S2 = sum(CYB(I,II+1:IEND)*SCR(II+1:IEND))
      S = S2+CYB(I,II)*SCR(II)
      S2 = CYB(I,II)+X*S2
      CYB(I,II) = S
      CYB(I,II+1:IEND) = CYB(I,II+1:IEND)-S2*SCR(II+1:IEND)
    end do
    X = One/CXA(II,II)
    do I=II+1,IEND
      CXA(I,II) = X*CXA(I,II)
      CXA(I,II+1:NTOT) = CXA(I,II+1:NTOT)-CXA(I,II)*CXA(II,II+1:NTOT)
    end do
    X = One/CYB(II,II)
    do I=II+1,IEND
      CYB(I,II) = X*CYB(I,II)
      CYB(I,II+1:NTOT) = CYB(I,II+1:NTOT)-CYB(I,II)*CYB(II,II+1:NTOT)
    end do
  end do
end do

end subroutine LU2
