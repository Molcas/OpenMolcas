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

subroutine Compute_M(ZA,nAtoms,RA,T,M)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtoms
real(kind=wp), intent(in) :: ZA(nAtoms), RA(3,nAtoms), T(3)
real(kind=wp), intent(out) :: M(3,3)
integer(kind=iwp) :: i, iAtom, j
real(kind=wp) :: RTx, RTy, RTz
real(kind=wp), parameter :: Thrs = 1.0e-14_wp

!                                                                      *
!***********************************************************************
!                                                                      *
! Form the nuclear charge moment tensor

M(:,:) = Zero
do iAtom=1,nAtoms
  RTx = RA(1,iAtom)-T(1)
  RTy = RA(2,iAtom)-T(2)
  RTz = RA(3,iAtom)-T(3)
  M(1,1) = M(1,1)+ZA(iAtom)*(RTy**2+RTz**2)
  M(2,2) = M(2,2)+ZA(iAtom)*(RTx**2+RTz**2)
  M(3,3) = M(3,3)+ZA(iAtom)*(RTx**2+RTy**2)

  M(1,2) = M(1,2)+ZA(iAtom)*(-RTx*RTy)
  M(1,3) = M(1,3)+ZA(iAtom)*(-RTx*RTz)
  M(2,1) = M(2,1)+ZA(iAtom)*(-RTy*RTx)

  M(2,3) = M(2,3)+ZA(iAtom)*(-RTy*RTz)
  M(3,1) = M(3,1)+ZA(iAtom)*(-RTz*RTx)
  M(3,2) = M(3,2)+ZA(iAtom)*(-RTz*RTy)
end do

! Remove noise

do i=1,3
  do j=1,3
    if (abs(M(i,j)) < Thrs) M(i,j) = Zero
  end do
end do
!call RecPrt('Compute_M: M',' ',M,3,3)

return

end subroutine Compute_M
