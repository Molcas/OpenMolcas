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

subroutine FIND_ACTIVE_BLOCKS(LUIN,LBLK,BLK_A,SEGMNT)
! Find the active (nonvanishing blocks) on LUIN
! Non vanishing block is flagged by a 1.0 ( note : real)
! in BLK_A

use lucia_data, only: IDISK

implicit none
integer LUIN, LBLK
! Output
real*8 BLK_A(*)
! Scratch
real*8 SEGMNT(*)
integer LBL(1), IDUMMY(1)
integer IBLK, NBLK_A, KBLK, NO_ZEROING, IMZERO, NBLK, NTEST, IAMPACK

IDISK(LUIN) = 0

IBLK = 0
NBLK_A = 0
! Loop over blocks
1000 continue
IBLK = IBLK+1
if (LBLK > 0) then
  LBL(1) = LBLK
else if (LBLK == 0) then
  call IDAFILE(LUIN,2,LBL,1,IDISK(LUIN))
else if (LBLK < 0) then
  call IDAFILE(LUIN,2,LBL,1,IDISK(LUIN))
  call IDAFILE(LUIN,2,IDUMMY,1,IDISK(LUIN))
end if
if (LBL(1) >= 0) then
  if (LBLK >= 0) then
    KBLK = LBL(1)
  else
    KBLK = -1
  end if
  NO_ZEROING = 1
  call FRMDSC2(SEGMNT,LBL(1),KBLK,LUIN,IMZERO,IAMPACK,NO_ZEROING)
  if (IMZERO == 0) then
    NBLK_A = NBLK_A+1
    BLK_A(IBLK) = 1.0d0
  else
    BLK_A(IBLK) = 0.0d0
  end if
end if
if ((LBL(1) >= 0) .and. (LBLK <= 0)) goto 1000
NBLK = IBLK-1

NTEST = 0
if (NTEST >= 1) write(6,*) ' FIND_A.... Number of total and active Blocks',NBLK,NBLK_A
if (NTEST >= 100) then
  write(6,*) ' Active blocks'
  call WRTMAT(BLK_A,1,NBLK,1,NBLK)
end if

end subroutine FIND_ACTIVE_BLOCKS
