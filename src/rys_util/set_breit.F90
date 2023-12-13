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
! Copyright (C) 2023, Roland Lindh                                     *
!***********************************************************************

subroutine Set_Breit(n)

use Breit, only: nOrdOp, nComp
use Definitions, only: iwp

integer(kind=iwp), intent(in) :: n

nOrdOp = n
if (nOrdOp == 0) then
  nComp = 1
else
  nComp = 6
end if

end subroutine Set_Breit

