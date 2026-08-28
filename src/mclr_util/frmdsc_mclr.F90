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

subroutine FRMDSC_MCLR(ARRAY,NDIM,MBLOCK,IFILE,IMZERO)
! TRANSFER ARRAY FROM DISC FILE IFILE

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NDIM, MBLOCK, IFILE
real(kind=wp), intent(out) :: ARRAY(NDIM)
integer(kind=iwp), intent(out) :: IMZERO
integer(kind=iwp) :: I, IBASE, IDUM(1), IREST, NBLOCK
integer(kind=iwp), parameter :: IPACK = 1, ICRAY = 1

if (IPACK /= 0) then
  ! Read if ARRAY is zero
  call IFRMDS(IDUM,1,MBLOCK,IFILE)
  IMZERO = IDUM(1)
  if (IMZERO == 1) then
    ARRAY(:) = Zero
    return
  end if
end if

if ((MBLOCK >= 0) .or. (ICRAY == 1)) then
  NBLOCK = MBLOCK
  if (MBLOCK <= 0) NBLOCK = NDIM
  IREST = NDIM
  IBASE = 0
  do while (IREST > 0)
    if (IREST > NBLOCK) then
      read(IFILE) (ARRAY(IBASE+I),I=1,NBLOCK)
      IBASE = IBASE+NBLOCK
      IREST = IREST-NBLOCK
    else
      read(IFILE) (ARRAY(IBASE+I),I=1,IREST)
      IREST = 0
    end if
  end do
end if

if ((MBLOCK < 0) .and. (NDIM > 0) .and. (ICRAY == 0)) then
  !call SQFILE(IFILE,2,ARRAY,2*NDIM)
  call SysHalt('frmdsc')
end if

return

end subroutine FRMDSC_MCLR
