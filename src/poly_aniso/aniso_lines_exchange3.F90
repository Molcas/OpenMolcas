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

subroutine Aniso_Lines_Exchange3(Jex,N1,N2,S1,S2,HAM)
! this Subroutine calculates the Lines exchange interaction between
! two sites, of the one interacting pair

use Constants, only: Zero, cZero, cOne
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: Jex(3)
integer(kind=iwp), intent(in) :: N1, N2
complex(kind=wp), intent(in) :: S1(3,N1,N1), S2(3,N2,N2)
complex(kind=wp), intent(out) :: HAM(N1,N1,N2,N2)
integer(kind=iwp) :: i1, i2, j1, j2, l
complex(kind=wp) :: Jc(3)
real(kind=wp), external :: dnrm2_

if ((N1 <= 0) .or. (N2 <= 0)) return
call zcopy_(N1*N1*N2*N2,[cZero],0,HAM,1)
if (dnrm2_(3,Jex,1) == Zero) return

call zcopy_(3,[cZero],0,Jc,1)
do l=1,3
  Jc(l) = -Jex(l)*cOne
end do

do i1=1,N1
  do j1=1,N1
    do i2=1,N2
      do j2=1,N2

        do l=1,3
          HAM(i1,j1,i2,j2) = HAM(i1,j1,i2,j2)+Jc(l)*S1(l,i1,j1)*S2(l,i2,j2)
        end do

      end do
    end do
  end do
end do

return

end subroutine Aniso_Lines_Exchange3
