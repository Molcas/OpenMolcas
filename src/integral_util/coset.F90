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

subroutine CoSet(iCoSet,nCoSet,iChAtom)

use Symmetry_Info, only: nIrrep, iOper
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: iCoSet(0:7), nCoSet
integer(kind=iwp), intent(in) :: iChAtom
integer(kind=iwp) :: iIrrep, itest, jCoSet, jtest
logical(kind=iwp) :: Same

! Find the coset representatives

iCoSet(0) = 0 ! Put in the unit operator
nCoSet = 1
do iIrrep=1,nIrrep-1
  itest = iand(iChAtom,iOper(iIrrep))
  Same = .false.
  do jCoSet=0,nCoSet-1
    jTest = iand(iChAtom,iCoSet(jCoSet))
    if (jTest == iTest) then
      Same = .true.
      exit
    end if
  end do
  if (.not. Same) then
    nCoSet = nCoSet+1
    iCoSet(nCoSet-1) = iOper(iIrrep)
  end if
end do

return

end subroutine CoSet
