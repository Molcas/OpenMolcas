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
  use strings
  implicit none

! This support up to about 100 active orbitals
  character(len=500), intent(inout) :: str
  character(len=1), dimension(500)  :: args
  character(len=500)                :: TempStr
  character(len=1) :: delims = ' '
  integer :: nargs, k, idxOrb

! Using George Benthien string utilites
  call parse(str, delims, args, nargs)

  idxOrb = 0

! Convert HFOC to Dice format
  str = ' '
  do k=1, nargs
    if (args(k) == '2') then
        write(TempStr, '(2I4)') idxOrb, idxOrb+1
        str = adjustl(trim(str)) // ' ' //adjustl(trim(TempStr))
    endif
    if (args(k) == 'u') then
        write(TempStr, '(I4)') idxOrb
        str = adjustl(trim(str)) // ' ' //adjustl(trim(TempStr))
    endif
    if (args(k) == 'd') then
        write(TempStr, '(I4)') idxOrb+1
        str = adjustl(trim(str)) // ' ' //adjustl(trim(TempStr))
    endif
    idxOrb = idxOrb + 2
  enddo

end subroutine
