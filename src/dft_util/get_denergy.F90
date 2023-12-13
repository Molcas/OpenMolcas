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

subroutine Get_dEnergy(Energy)

use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: Energy
logical(kind=iwp) :: Found_EAV

Found_EAV = .false.
call Qpg_dScalar('Average energy',Found_EAV)

if (Found_EAV) then
  call Get_dScalar('Average energy',Energy)
else
  call Get_dScalar('Last energy',Energy)
end if

return

end subroutine Get_dEnergy
