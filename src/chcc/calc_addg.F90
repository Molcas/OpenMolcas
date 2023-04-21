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

subroutine Calc_addG(aGrp,adda)
! calc add constant from group

use chcc_global, only: DimGrpa
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: aGrp
integer(kind=iwp), intent(out) :: adda

adda = sum(DimGrpa(1:aGrp-1))

return

end subroutine Calc_addG
