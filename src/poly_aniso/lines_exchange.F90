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

subroutine Lines_Exchange(Jex,N1,N2,S1,S2,HAM)
! this Subroutine calculates the Lines exchange interaction between
! two sites, of the one interacting pair

use Constants, only: Zero, cZero

implicit none
! input variables
integer, intent(in) :: N1, N2
real(kind=8), intent(in) :: Jex
complex(kind=8), intent(in) :: S1(3,N1,N1)
complex(kind=8), intent(in) :: S2(3,N2,N2)
! output variables
complex(kind=8), intent(out) :: HAM(N1,N1,N2,N2)
! local variables
integer :: i1, i2, j1, j2, l

if ((N1 <= 0) .or. (N2 <= 0)) return
call zcopy_(N1*N1*N2*N2,[cZero],0,HAM,1)
if (Jex == Zero) return

! kind=8, complex double precision
do i1=1,N1
  do j1=1,N1
    do i2=1,N2
      do j2=1,N2

        do l=1,3
          HAM(i1,j1,i2,j2) = HAM(i1,j1,i2,j2)-Jex*S1(l,i1,j1)*S2(l,i2,j2)
        end do

      end do
    end do
  end do
end do

return

end subroutine Lines_Exchange
