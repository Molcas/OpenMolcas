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
! Copyright (C) 1990, IBM                                              *
!***********************************************************************

function NrOpr(iOp)
!***********************************************************************
!                                                                      *
! Object: to return the order index of a symmetry operator.            *
!                                                                      *
!***********************************************************************

use Symmetry_Info, only: iOper, nIrrep
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: NrOpr
integer(kind=iwp), intent(in) :: iOp
integer(kind=iwp) :: iIrrep

NrOpr = -1
do iIrrep=0,nIrrep-1
  if (iOp == iOper(iIrrep)) NrOpr = iIrrep
end do

return

end function NrOpr
