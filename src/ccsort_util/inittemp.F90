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

subroutine inittemp(num)
! this routine initializes status matrix
! num - number of files to be used (I)

use ccsort_global, only: lrectemp, nrectemp, stattemp
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: num

stattemp(1:num) = 0
nrectemp(1:num) = 0
lrectemp(1:num) = 0

return

end subroutine inittemp
