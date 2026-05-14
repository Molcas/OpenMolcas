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

subroutine SKPRCD2(NDIM,MBLOCK,IFILE)
! Skip record in file IFILE
!
! Version allowing zero and packed blocks
!
! Does not work with FASTIO - I expect

use lucia_data, only: IDISK
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NDIM, MBLOCK, IFILE
integer(kind=iwp) :: I_AM_PACKED, IBASE, IDUMMY(1), IMZERO, IPACK, IREST, ISCR(2), ISTOP, LBATCH, NBLOCK
real(kind=wp) :: DUMMY(1)

IPACK = 1
if (IPACK /= 0) then
  ! Read if ARRAY is zero
  call IFRMDS(ISCR,2,2,IFILE)
  IMZERO = ISCR(1)
  I_AM_PACKED = ISCR(2)
  if (IMZERO == 1) return
end if

if (I_AM_PACKED == 1) then
  ! Loop over packed records of dimension LPBLK
  ! The next LPBLK elements
  ISTOP = 0
  do while (ISTOP == 0)
    ! Read next batch
    call IDAFILE(IFILE,2,ISCR,1,IDISK(IFILE))
    LBATCH = ISCR(1)
    if (LBATCH > 0) then
      IDUMMY(1) = 0
      call IDAFILE(IFILE,0,IDUMMY,LBATCH,IDISK(IFILE))
      DUMMY(1) = Zero
      call DDAFILE(IFILE,0,DUMMY,LBATCH,IDISK(IFILE))
    end if
    call IDAFILE(IFILE,2,ISCR,1,IDISK(IFILE))
    ISTOP = ISCR(1)
  end do
else if (I_AM_PACKED == 0) then
  !vv if (.true.) then
  NBLOCK = MBLOCK
  if (MBLOCK <= 0) NBLOCK = NDIM
  IREST = NDIM
  IBASE = 0
  do
    if (IREST > NBLOCK) then
      call DDAFILE(IFILE,0,DUMMY,NBLOCK,IDISK(IFILE))
      IBASE = IBASE+NBLOCK
      IREST = IREST-NBLOCK
    else
      call DDAFILE(IFILE,0,DUMMY,IREST,IDISK(IFILE))
      IREST = 0
    end if
    call IDAFILE(IFILE,0,IDUMMY,1,IDISK(IFILE))
    if (IREST <= 0) exit
  end do
  !vv end if

end if

end subroutine SKPRCD2
