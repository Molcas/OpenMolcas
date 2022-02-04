!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Valera Veryazov                                        *
!***********************************************************************
!  isFreeUnit
!
!> @brief
!>   Find free unit number
!> @author V. Veryazov
!>
!> @details
!> Find unused unit number, starting from initial value.
!>
!> @param[in] iseed guess for unit number
!>
!> @return Free unit number
!***********************************************************************

function isFreeUnit(iseed)
! check free chanel, starting from iseed

use Fast_IO, only: isOpen, MxFile
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: isFreeUnit
integer(kind=iwp), intent(in) :: iseed
integer(kind=iwp) :: init, kan, kan0
logical(kind=iwp) :: is_opened, skip

!VV since more and more developers' calling isfreeunit with constant...
init = iseed
if ((init < 1) .or. (init > 300)) then
  write(u6,*) '*** Possible bug in opening file'
  write(u6,*) '*** isFreeUnit resets the unit number'
  write(u6,*) 'init=',init
  init = 12
  call Abend()
end if
isFreeUnit = -init
kan = min(init,MxFile-1)
kan0 = kan
skip = .true.
do
  if (.not. skip) then
    if (kan == MxFile+1) kan = 10
    if (kan == kan0) then
      call fastio('STATUS')
      write(u6,*) ' isFreeUnit: no available unit!'
      call Abend()
    end if
  end if
  skip = .false.

  ! Check for Dafile

  if ((kan > 1) .and. (kan <= MxFile) .and. (isOpen(kan) == 1)) then
    kan = kan+1
  else

    ! Check for Fortran I/O

    inquire(unit=kan,opened=is_opened)
    if (is_opened) then
      kan = kan+1
    else
      exit
    end if
  end if
end do
isFreeUnit = kan

return

end function isFreeUnit
