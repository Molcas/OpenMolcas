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

function nU(iU)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: nU
integer(kind=iwp), intent(in) :: iU
integer(kind=iwp) :: i

nU = 0
do i=0,7
  if (btest(iU,i)) nU = nU+1
end do

return

end function nU
