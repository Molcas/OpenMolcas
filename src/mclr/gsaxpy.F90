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
if (IWAY == 1) then
  ! Straightforward sequence
  do J=1,NABCOL
    do K=1,NACOL
      JACT = IABCOL(J)
      KACT = IACOL(K)
      AB(:,JACT) = AB(:,JACT)+B(K,J)*A(:,KACT)
    end do
  end do

else if (IWAY == 2) then
  ! Unrolling over columns of A
  NROL = 5
  NRES = mod(NACOL,NROL)
  ! overhead
  select case (NRES)
    case (1)
      do J=1,NABCOL
        K1ACT = IACOL(1)
        JACT = IABCOL(J)
        AB(:,JACT) = AB(:,JACT)+B(1,J)*A(:,K1ACT)
      end do
    case (2)
      do J=1,NABCOL
        K1ACT = IACOL(1)
        K2ACT = IACOL(2)
        JACT = IABCOL(J)
        AB(:,JACT) = AB(:,JACT)+B(1,J)*A(:,K1ACT)+B(2,J)*A(:,K2ACT)
      end do
    case (3)
      do J=1,NABCOL
        K1ACT = IACOL(1)
        K2ACT = IACOL(2)
        K3ACT = IACOL(3)
        JACT = IABCOL(J)
        AB(:,JACT) = AB(:,JACT)+B(1,J)*A(:,K1ACT)+B(2,J)*A(:,K2ACT)+B(3,J)*A(:,K3ACT)
      end do
    case (4)
      do J=1,NABCOL
        K1ACT = IACOL(1)
        K2ACT = IACOL(2)
        K3ACT = IACOL(3)
        K4ACT = IACOL(4)
        JACT = IABCOL(J)
        AB(:,JACT) = AB(:,JACT)+B(1,J)*A(:,K1ACT)+B(2,J)*A(:,K2ACT)+B(3,J)*A(:,K3ACT)+B(4,J)*A(:,K4ACT)
      end do
    case default
  end select
  ! (End of Overhead)
  do K=NRES+1,NACOL,NROL
    do J=1,NABCOL
      K1ACT = IACOL(K)
      K2ACT = IACOL(K+1)
      K3ACT = IACOL(K+2)
      K4ACT = IACOL(K+3)
      K5ACT = IACOL(K+4)
      JACT = IABCOL(J)
      AB(:,JACT) = AB(:,JACT)+B(K,J)*A(:,K1ACT)+B(K+1,J)*A(:,K2ACT)+B(K+2,J)*A(:,K3ACT)+B(K+3,J)*A(:,K4ACT)+B(K+4,J)*A(:,K5ACT)
    end do
  end do
end if
! (End of IWAY branching)

return

end subroutine GSAXPY
