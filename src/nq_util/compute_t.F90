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

subroutine Compute_T(Z_Tot,T,ZA,RA,nAtoms)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtoms
real(kind=wp), intent(in) :: Z_Tot, ZA(nAtoms), RA(3,nAtoms)
real(kind=wp), intent(out) :: T(3)
integer(kind=iwp) :: iAtom, iCar
real(kind=wp) :: Tmp

!                                                                      *
!***********************************************************************
!                                                                      *
! Form the center of nuclear charge

do iCar=1,3
  Tmp = Zero
  do iAtom=1,nAtoms
    Tmp = Tmp+ZA(iAtom)*RA(iCar,iAtom)
  end do
  T(iCar) = Tmp/Z_Tot
end do
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Compute_T
