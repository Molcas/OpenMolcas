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

function CovRadT(i)

use CovRad_Data, only: CovRadT_
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: CovRadT
integer(kind=iwp), intent(in) :: i

if (i > 92) then
  !write(u6,*) 'CovRadT: i > 92'
  CovRadT = 1.50_wp
else
  CovRadT = CovRadT_(i)
end if

return

end function CovRadT
