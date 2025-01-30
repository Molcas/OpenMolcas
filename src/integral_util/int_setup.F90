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

subroutine Int_Setup(Coor)

use Gateway_Info, only: DoFMM, RPQMin
use Gateway_global, only: FMM_shortrange
use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: Coor(3,4)
integer(kind=iwp) :: i
real(kind=wp) :: D, P, Q

! For the FMM coulomb integrals <AB(r1)|1/r12|CD(r2)>
! Here we flag the integral routines that we only want to compute
! the short-range non-multipole component of integrals over this
! shell quartet if midpoint(A,B) is sufficiently far from
! midpoint(C,D) for numerical stability.
! Note that midpoint(X,Y) corresponds to the multipole expansion
! centre of an XY AO-pair, regardless of exponents.

FMM_shortrange = .false.
if (DoFMM) then
  D = Zero
  do i=1,3
    P = Half*(Coor(i,1)+Coor(i,2))    ! AB shell-pair
    Q = Half*(Coor(i,3)+Coor(i,4))    ! CD shell-pair
    D = D+(P-Q)*(P-Q)
  end do
  if (D > RPQMIN*RPQMIN) FMM_shortrange = .true.
end if

end subroutine Int_Setup
