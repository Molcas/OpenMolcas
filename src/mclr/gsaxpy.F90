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
! Copyright (C) 1994, Jeppe Olsen                                      *
!***********************************************************************

subroutine GSAXPY(AB,A,B,NABCOL,NACOL,NROW,IABCOL,IACOL)
! AB(I,IABCOL(J)) = AB(I,IABCOL(J)) + A(I,IACOL(K))*B(K,J)
!
! Jeppe Olsen, Spring of 94 Daughter of MSAXPY*

implicit real*8(A-H,O-Z)
! Input
dimension A(NROW,*), B(NACOL,NABCOL)
dimension IACOL(NACOL), IABCOL(NABCOL)
! Output
dimension AB(NROW,*)

IWAY = 2
!ICRAY = 0
if (IWAY == 1) then
  ! Straightforward sequence of SAXPY's
  do J=1,NABCOL
    do K=1,NACOL
      JACT = IABCOL(J)
      KACT = IACOL(K)
      FACTOR = B(K,J)
      !if (ICRAY == 1) then
      !  call SAXPY(NROW,FACTOR,A(1,KACT),1,AB(1,JACT),1)
      !else
      do I=1,NROW
        AB(I,JACT) = AB(I,JACT)+FACTOR*A(I,KACT)
      end do
      !end if
    end do
  end do

else if (IWAY == 2) then
  ! Unrolling over columns of A
  NROL = 5
  NRES = mod(NACOL,NROL)
  ! overhead
  if (NRES == 1) then
    do J=1,NABCOL
      JACT = IABCOL(J)
      K1ACT = IACOL(1)
      B1J = B(1,J)
      !if (ICRAY == 1) then
      !  call SAXPY(NROW,B1J,A(1,K1ACT),1,AB(1,JACT),1)
      !else
      do I=1,NROW
        AB(I,JACT) = AB(I,JACT)+A(I,K1ACT)*B1J
      end do
      !endif
    end do
  else if (NRES == 2) then
    do J=1,NABCOL
      K1ACT = IACOL(1)
      K2ACT = IACOL(2)
      B1J = B(1,J)
      B2J = B(2,J)
      JACT = IABCOL(J)
      do I=1,NROW
        AB(I,JACT) = AB(I,JACT)+A(I,K1ACT)*B1J+A(I,K2ACT)*B2J
      end do
    end do
  else if (NRES == 3) then
    do J=1,NABCOL
      K1ACT = IACOL(1)
      K2ACT = IACOL(2)
      K3ACT = IACOL(3)
      JACT = IABCOL(J)
      B1J = B(1,J)
      B2J = B(2,J)
      B3J = B(3,J)
      do I=1,NROW
        AB(I,JACT) = AB(I,JACT)+A(I,K1ACT)*B1J+A(I,K2ACT)*B2J+A(I,K3ACT)*B3J
      end do
    end do
  else if (NRES == 4) then
    do J=1,NABCOL
      K1ACT = IACOL(1)
      K2ACT = IACOL(2)
      K3ACT = IACOL(3)
      K4ACT = IACOL(4)
      JACT = IABCOL(J)
      B1J = B(1,J)
      B2J = B(2,J)
      B3J = B(3,J)
      B4J = B(4,J)
      do I=1,NROW
        AB(I,JACT) = AB(I,JACT)+A(I,K1ACT)*B1J+A(I,K2ACT)*B2J+A(I,K3ACT)*B3J+A(I,K4ACT)*B4J
      end do
    end do
  end if
  ! (End of Overhead)
  do K=NRES+1,NACOL,NROL
    do J=1,NABCOL
      K1ACT = IACOL(K)
      K2ACT = IACOL(K+1)
      K3ACT = IACOL(K+2)
      K4ACT = IACOL(K+3)
      K5ACT = IACOL(K+4)
      JACT = IABCOL(J)
      B1J = B(K,J)
      B2J = B(K+1,J)
      B3J = B(K+2,J)
      B4J = B(K+3,J)
      B5J = B(K+4,J)
      do I=1,NROW
        AB(I,JACT) = AB(I,JACT)+A(I,K1ACT)*B1J+A(I,K2ACT)*B2J+A(I,K3ACT)*B3J+A(I,K4ACT)*B4J+A(I,K5ACT)*B5J
      end do
    end do
  end do
end if
! (End of IWAY branching)

return

end subroutine GSAXPY
