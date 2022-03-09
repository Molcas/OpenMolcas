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

subroutine FMMM(A,B,C,NROW,NCOL,NSUM)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NROW, NCOL, NSUM
real(kind=wp), intent(in) :: A(NROW,NSUM), B(NSUM,NCOL)
real(kind=wp), intent(out) :: C(NROW,NCOL)
integer(kind=iwp) :: I, J, K, KK
real(kind=wp) :: T
integer(kind=iwp), parameter :: KS = 48

C(:,:) = Zero

do KK=1,NSUM,KS
  do I=1,NROW
    do J=1,NCOL
      T = C(I,J)
      do K=KK,min(KK+KS-1,NSUM)
        T = T+B(K,J)*A(I,K)
      end do
      C(I,J) = T
    end do
  end do
end do

return

end subroutine FMMM
