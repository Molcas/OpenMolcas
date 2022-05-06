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

subroutine Truncate_Grid(R,nR,nR_Eff,Radius_Max)

implicit real*8(a-h,o-z)
real*8 R(2,nR)

nTmp = nR_Eff
do i=1,nTmp
  if (R(1,i) > Radius_Max) then
    nR_Eff = i-1
    Go To 99
  end if
end do
99 continue

return

end subroutine Truncate_Grid
