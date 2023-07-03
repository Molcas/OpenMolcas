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

subroutine updRotMat(SmallRot,RotMat)

implicit real*8(a-h,o-z)
dimension SmallRot(3), RotMat(3,3)
dimension tmp(3,3), SmallMat(3,3)

call mkRotMat(SmallRot,SmallMat)
do i=1,3
  do j=1,3
    sum = 0.0d0
    do k=1,3
      sum = sum+RotMat(i,k)*SmallMat(k,j)
    end do
    tmp(i,j) = sum
  end do
end do
do i=1,3
  do j=1,3
    RotMat(i,j) = tmp(i,j)
  end do
end do
! Check for orthonormality:
do i=1,3
  do j=1,3
    sum = 0.0d0
    if (i == j) sum = -1.0d0
    do k=1,3
      sum = sum+RotMat(i,k)*RotMat(j,k)
    end do
    if (abs(sum) > 1.0D-10) then
      write(6,*) ' UPDROTMAT ON check sum error:',sum
    end if
  end do
end do

return

end subroutine updRotMat
