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

subroutine driverBas(ireturn)

use desymmetrize_mod, only: desym
use info_expbas_mod, only: DoDesy, DoExpbas
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: ireturn

ireturn = 0

call Readinp_expbas()
if (DoExpbas) then
  call expbas(ireturn)
  if (ireturn /= 0) return
end if

if (DoDesy) then
  call desym(ireturn)
  if (ireturn /= 0) return
end if

end subroutine driverBas
