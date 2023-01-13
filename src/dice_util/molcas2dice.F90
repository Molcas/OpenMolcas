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
! Copyright (C) 2022, Quan Phung                                       *
!***********************************************************************
! Subroutine to convert HFOC to Dice format
! Written by Quan Phung, Leuven, Apr 2018
!                        Nagoya, Oct 2022

subroutine molcas2dice(str)

use Definitions, only: iwp

implicit none
! This supports up to about 100 active orbitals
character(len=500), intent(inout) :: str
character(len=500) :: str2
character(len=8) :: TempStr
integer(kind=iwp) :: i, idxOrb

idxOrb = 0

! Convert DIOC to Dice format
str2 = ' '
do i=1,len(str)
  select case (str(i:i))
    case ('2')
      write(TempStr,'(2I4)') idxOrb,idxOrb+1
      str2 = adjustl(trim(str2))//' '//adjustl(trim(TempStr))
      idxOrb = idxOrb+2
    case ('u','a')
      write(TempStr,'(I4)') idxOrb
      str2 = adjustl(trim(str2))//' '//adjustl(trim(TempStr))
      idxOrb = idxOrb+2
    case ('d','b')
      write(TempStr,'(I4)') idxOrb
      str2 = adjustl(trim(str2))//' '//adjustl(trim(TempStr))
      idxOrb = idxOrb+2
    case (' ')
      ! do nothing
    case default
      idxOrb = idxOrb+2
  end select
end do

str = str2

end subroutine molcas2dice
