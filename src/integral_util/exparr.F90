!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

subroutine ExpArr(Array,Ind,nArray,lArray)
!***********************************************************************
!                                                                      *
! Object: to expand arrays according to an index array.                *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             August '91                                               *
!***********************************************************************

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nArray, Ind(nArray), lArray
real(kind=wp), intent(inout) :: Array(lArray,nArray)
integer(kind=iwp) :: iArray, jArray

do iArray=nArray,1,-1
  jArray = Ind(iArray)
  if (jArray <= 0) then
    ! Set column iArray to zero
    Array(:,iArray) = Zero
  else if (jArray < iArray) then
    ! Copy row jArray to position iArray
    Array(:,iArray) = Array(:,jArray)
  end if
end do

return

end subroutine ExpArr
