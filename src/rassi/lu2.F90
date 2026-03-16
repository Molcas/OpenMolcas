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
!  GIVES A SIMULTANEOUS LU-PARTITIONING OF MATRICES CXA,CYB IN THE
!  SENSE THAT CXA*X = L1*U1 AND CYB*X = L2*U2, WHERE X IS A BLOCK
!  UNITARY MATRIX. X IS NEVER FORMED, BUT AT EACH STEP SUCH A
!  TRANSFORMATION IS APPLIED THAT THE LU PARTITIONING CAN BE
!  APPLIED IN AN OPTIMAL WAY -- THIS IS REFERRED TO AS A
!  UNITARY PSEUDO-PIVOTATION. MATRICES CXA AND CYB ARE DESTROYED,
!  AND WILL CONTAIN THE NONTRIVIAL ELEMENTS OF THE TRIANGULAR
!  MATRICES.
!                                         ( MALMQUIST 84-01-16 )

use definitions, only: iwp, wp, u6
use constants, only: Zero, One, Two

implicit none
integer(kind=iwp), intent(in) :: NDIMEN, NBLOCK
integer(kind=iwp), intent(in) :: NSIZE(NBLOCK)
real(kind=wp), intent(inout) :: CXA(NDIMEN,NDIMEN), CYB(NDIMEN,NDIMEN)
real(kind=wp), intent(inout) :: SCR(NDIMEN)
integer(kind=iwp) NTOT, IBLOCK, IEND, ISTA, II, J, I
real(kind=wp) S1, S2, S3, X1, X2, S, X
real(kind=wp), parameter :: THR = 1.0e-6_wp

NTOT = 0
do IBLOCK=1,NBLOCK
  NTOT = NTOT+NSIZE(IBLOCK)
end do
if (NTOT /= NDIMEN) write(u6,*) ' ERROR: NTOT /= NDIMEN IN LU2.'
IEND = 0
do IBLOCK=1,NBLOCK
  ISTA = IEND+1
  IEND = IEND+NSIZE(IBLOCK)
  do II=ISTA,IEND
    S1 = Zero
    S2 = Zero
    S3 = Zero
    do J=II,IEND
      S1 = S1+CXA(II,J)**2
      S2 = S2+CYB(II,J)**2
      S3 = S3+CXA(II,J)*CYB(II,J)
    end do
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
    do I=II,IEND
      SCR(I) = X1*CXA(II,I)+X2*CYB(II,I)
    end do
    S = Two*(One+X1*X2*S3)
    X = One/sign(sqrt(S),SCR(II))
    do I=II,IEND
      SCR(I) = X*SCR(I)
    end do
    X = One/(One+SCR(II))
    do I=1,IEND
      S = Zero
      do J=II,IEND
        S = S+CXA(I,J)*SCR(J)
      end do
      S2 = Zero
      do J=II+1,IEND
        S2 = S2+CXA(I,J)*SCR(J)
      end do
      S2 = CXA(I,II)+X*S2
      CXA(I,II) = S
      do J=II+1,IEND
        CXA(I,J) = CXA(I,J)-S2*SCR(J)
      end do
    end do
    do I=1,IEND
      S = Zero
      do J=II,IEND
        S = S+CYB(I,J)*SCR(J)
      end do
      S2 = Zero
      do J=II+1,IEND
        S2 = S2+CYB(I,J)*SCR(J)
      end do
      S2 = CYB(I,II)+X*S2
      CYB(I,II) = S
      do J=II+1,IEND
        CYB(I,J) = CYB(I,J)-S2*SCR(J)
      end do
    end do
    X = One/CXA(II,II)
    do I=II+1,IEND
      CXA(I,II) = X*CXA(I,II)
      do J=II+1,NTOT
        CXA(I,J) = CXA(I,J)-CXA(I,II)*CXA(II,J)
      end do
    end do
    X = One/CYB(II,II)
    do I=II+1,IEND
      CYB(I,II) = X*CYB(I,II)
      do J=II+1,NTOT
        CYB(I,J) = CYB(I,J)-CYB(I,II)*CYB(II,J)
      end do
    end do
  end do
end do

end subroutine LU2
