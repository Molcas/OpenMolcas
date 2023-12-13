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

subroutine Energy_E2d(V,Tau,e,eos,dima,no)
! this routine calcs:
! E2 = sum(a,b,i,j) (2V(ai|bj)-V(aj|bi)) . Tau(a,b,i,j
! E2os=sum(a,b,i,j) (  (ai|bj)         ) . Tau(a,b,i,j
! where
!    E2    - complete E2 component of energy
!    E2os  - other spin E2 component of energy
! for aGrp=bGrp (Tau array is full, but only a>=b values are completed)

use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, no
real(kind=wp), intent(in) :: V(dima,no,dima,no), Tau(dima,dima,no,no)
real(kind=wp), intent(out) :: e, eos
integer(kind=iwp) :: a, b, i, j
real(kind=wp) :: ehlp

e = Zero
eos = Zero
ehlp = Zero

! off diagonal

do j=1,no
  do i=1,no
    do b=1,dima-1
      ehlp = ehlp+V(b,i,b,j)*Tau(b,b,i,j)
      do a=b+1,dima
        e = e+(Two*V(a,i,b,j)-V(a,j,b,i))*Tau(a,b,i,j)
        eos = eos+V(a,i,b,j)*Tau(a,b,i,j)
      end do
    end do
    ehlp = ehlp+V(dima,i,dima,j)*Tau(dima,dima,i,j)
  end do
end do

e = Two*e+ehlp
eos = Two*eos+ehlp

return

end subroutine Energy_E2d
