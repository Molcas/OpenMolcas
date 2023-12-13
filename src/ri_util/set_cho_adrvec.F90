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

function Set_CHO_ADRVEC(ii)

use Cholesky, only: CHO_ADRVEC
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: Set_CHO_ADRVEC
integer(kind=iwp), intent(in) :: ii

Set_CHO_ADRVEC = 0
if (ii < 0) then
  Set_CHO_ADRVEC = CHO_ADRVEC
else if ((ii == 1) .or. (ii == 2)) then
  CHO_ADRVEC = ii
  Set_CHO_ADRVEC = CHO_ADRVEC
else
  call WarningMessage(2,'Set_CHO_ADRVEC: Illegal option')
  write(u6,*) 'ii=',ii
  call Abend()
end if

return

end function Set_CHO_ADRVEC
