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

subroutine Get_bastype(BasisTypes,nData)

implicit real*8(A-H,O-Z)
integer BasisTypes(4), loc_BasisTypes(4)
integer is_BasisType
save is_BasisType
data is_BasisType/0/
save loc_BasisTypes

!Label = 'BasType'
if (is_BasisType == 0) then
  call get_iArray('BasType',loc_BasisTypes,nData)
  is_BasisType = 1
end if
do i=1,nData
  BasisTypes(i) = loc_BasisTypes(i)
end do

return

end subroutine Get_bastype
