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

subroutine CalcAddpp(aSGrp,addapp)
! this routine calcs addapp

use chcc_global, only: DimSGrpa
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: aSGrp
integer(kind=iwp), intent(out) :: addapp

addapp = 0
if (aSGrp > 1) addapp = sum(DimSGrpa(1:aSGrp-1))

return

end subroutine CalcAddpp
