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

subroutine IFRMDS(IARRAY,NDIM,MBLOCK,IFILE)
! TRANSFER INTEGER ARRAY FROM DISC FILE IFILE
!
! NBLOCK < 0 INDICATES USE OF FASTIO
!
! If nblock == 0 NBLOCK = NDIM

use lucia_data, only: IDISK
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: NDIM, MBLOCK, IFILE
integer(kind=iwp), intent(out) :: IARRAY(NDIM)
integer(kind=iwp) :: IBASE, IDUMMY(1), IREST, NBLOCK

NBLOCK = MBLOCK

! DO NOT USE FASTIO
if (NBLOCK <= 0) NBLOCK = NDIM
IREST = NDIM
IBASE = 0
do
  if (IREST > NBLOCK) then
    call IDAFILE(IFILE,2,IARRAY(IBASE+1),NBLOCK,IDISK(IFILE))
    IBASE = IBASE+NBLOCK
    IREST = IREST-NBLOCK
  else
    call IDAFILE(IFILE,2,IARRAY(IBASE+1),IREST,IDISK(IFILE))
    IREST = 0
  end if
  call IDAFILE(IFILE,2,IDUMMY,1,IDISK(IFILE))
  if (IREST <= 0) exit
end do

end subroutine IFRMDS
