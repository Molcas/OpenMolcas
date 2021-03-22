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

function RPA_iUHF()

use RPA_globals, only: Reference
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: RPA_iUHF
integer(kind=iwp) :: iUHF

if (Reference(1:1) == 'R') then
  iUHF = 1
else if (Reference(1:1) == 'U') then
  iUHF = 2
else
  write(u6,'(A,A)') 'Reference=',Reference
  call RPA_Warn(3,'Unable to determine iUHF in RPA')
  iUHF = -1
end if
RPA_iUHF = iUHF

end function RPA_iUHF
