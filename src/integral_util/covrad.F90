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

function CovRad(i)

use CovRad_Data, only: CovRad_
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: CovRad
integer(kind=iwp), intent(in) :: i

if (i > 86) then
  !write(Warning,'(2A)') 'CovRad: Warning i > 86!,;Guesstimate of 2.70 au is used!'
  !call WarningMessage(1,Warning)
  CovRad = 2.70_wp
else
  CovRad = CovRad_(i)
end if

return

end function CovRad
