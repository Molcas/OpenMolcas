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
! NOTE: THE FOLLOWING ROUTINE MAY NOT BE VECTORIZED.
! THERE IS A COMPILER BUG ON FORTRAN VERSION 2.2.0 (JUNE 1987).
! WRONG RESULTS PRODUCED EVEN ON VECTORIZED LEVEL 1.

subroutine MKCXAL(NDIMEN,TRAL,CXAL)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NDIMEN
real(kind=wp), intent(in) :: TRAL(NDIMEN,NDIMEN)
real(kind=wp), intent(out) :: CXAL(NDIMEN,NDIMEN)
integer(kind=iwp) :: I, J, K
real(kind=wp) :: SUMMA

do I=1,NDIMEN
  do J=I,NDIMEN
    CXAL(I,J) = Zero
  end do
  CXAL(I,I) = One
end do
do K=1,NDIMEN
  do I=1,K-1
    SUMMA = Zero
    do J=1,K-1
      SUMMA = SUMMA+CXAL(I,J)*TRAL(J,K)
    end do
    CXAL(I,K) = -(SUMMA/TRAL(K,K))
  end do
  do I=K,NDIMEN
    SUMMA = TRAL(I,K)
    if (I == K) SUMMA = -One
    do J=1,K-1
      SUMMA = SUMMA+CXAL(I,J)*TRAL(J,K)
    end do
    CXAL(I,K) = -(SUMMA/TRAL(K,K))
  end do
end do

end subroutine MKCXAL
