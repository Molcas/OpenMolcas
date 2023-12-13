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

subroutine Union(iU,nU,iV,nV,iR,iM,nM)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nU, iU(nU), nV, iV(nV), iR
integer(kind=iwp), intent(out) :: iM(8), nM
integer(kind=iwp) :: i, iRV
logical(kind=iwp) :: RinT_

! M is formed as U union RU

iM(1:nU) = iU(:) ! copy the first elements
nM = nU
do i=1,nV
  iRV = ieor(iR,iV(i))
  if (.not. RinT_(iM,nM,iRV)) then
    nM = nM+1
    iM(nM) = iRV
  end if
end do

return

end subroutine Union
