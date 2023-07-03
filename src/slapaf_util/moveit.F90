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

implicit real*8(a-h,o-z)
real*8 U(3,3)
real*8 CurrXYZ(3,nMass)
real*8 Ref123(3,nMass)
real*8 RotVec(3)
real*8 Trans(3)

call mkRotMat(RotVec,U)
do imass=1,nmass
  do i=1,3
    sum = Trans(i)
    do j=1,3
      sum = sum+U(i,j)*Ref123(j,imass)
    end do
    CurrXYZ(i,imass) = sum
  end do
end do

return

end subroutine moveit
