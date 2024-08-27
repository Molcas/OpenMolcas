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

implicit none
integer, intent(out) :: nCoSet
integer, intent(In) :: iChAtom
integer, intent(Out) :: iCoSet(0:7)
logical Same
integer iIrrep, itest, jCoSet, jtest

! Find the coset representatives

iCoSet(0) = 0 ! Put in the unit operator
nCoSet = 1
do iIrrep=1,nIrrep-1
  itest = iand(iChAtom,iOper(iIrrep))
  Same = .false.
  do jCoSet=0,nCoSet-1
    jTest = iand(iChAtom,iCoSet(jCoSet))
    Same = Same .or. (jTest == iTest)
  end do
  if (.not. Same) then
    nCoSet = nCoSet+1
    iCoSet(nCoSet-1) = iOper(iIrrep)
  end if
end do

return

end subroutine CoSet
