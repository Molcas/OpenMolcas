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

function iNew(iTest,iIrrep)

use Symmetry_Info, only: iChTbl, nIrrep
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iNew
integer(kind=iwp), intent(in) :: iIrrep, iTest(8)
integer(kind=iwp) :: i, iGo, j

iNew = 0
! Test iTest against all rows thus far.
do i=1,iIrrep
  ! Do a intger inner product.
  iGo = 0
  do j=1,nIrrep
    iGo = iGo+iTest(j)*iChTbl(i-1,j-1)
  end do
  if (iGo /= 0) then
    ! Here if row already defined.
    iNew = i
    return
  end if
end do
iNew = iIrrep+1

return

end function iNew
