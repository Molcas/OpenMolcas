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

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NDIMEN
real(kind=wp), intent(in) :: TRAL(NDIMEN,NDIMEN)
real(kind=wp), intent(out) :: CXAL(NDIMEN,NDIMEN)
integer(kind=iwp) :: I, K
real(kind=wp) :: SUMMA

call unitmat(CXAL,NDIMEN)
do K=1,NDIMEN
  do I=1,K-1
    SUMMA = sum(CXAL(I,1:K-1)*TRAL(1:K-1,K))
    CXAL(I,K) = -SUMMA/TRAL(K,K)
  end do
  do I=K,NDIMEN
    if (I == K) then
      SUMMA = -One
    else
      SUMMA = TRAL(I,K)
    end if
    SUMMA = SUMMA+sum(CXAL(I,1:K-1)*TRAL(1:K-1,K))
    CXAL(I,K) = -SUMMA/TRAL(K,K)
  end do
end do

end subroutine MKCXAL
