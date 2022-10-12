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

subroutine Qpg_zArray(Label,Found,nData)

use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: Label
logical(kind=iwp), intent(out) :: Found
integer(kind=iwp), intent(out) :: nData
integer(kind=iwp) :: nDataI, nDataR
logical(kind=iwp) :: FoundI, FoundR

call qpg_darray('R'//Label,FoundR,nDataR)
call qpg_darray('I'//Label,FoundI,nDataI)

if ((nDataR == nDataI) .and. FoundR .and. FoundI) then
  nData = nDataR
  found = .true.
else
  found = .false.
  nData = 0
end if

return

end subroutine Qpg_zArray
