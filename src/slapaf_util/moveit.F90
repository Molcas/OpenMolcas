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

subroutine moveit(nmass,currxyz,ref123,trans,rotvec)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nMass
real(kind=wp), intent(out) :: CurrXYZ(3,nMass)
real(kind=wp), intent(in) :: Ref123(3,nMass), Trans(3), RotVec(3)
real(kind=wp) :: U(3,3)
integer(kind=iwp) :: i, imass, j
real(kind=wp) :: rsum

call mkRotMat(RotVec,U)
do imass=1,nmass
  do i=1,3
    rsum = Trans(i)
    do j=1,3
      rsum = rsum+U(i,j)*Ref123(j,imass)
    end do
    CurrXYZ(i,imass) = rsum
  end do
end do

return

end subroutine moveit
