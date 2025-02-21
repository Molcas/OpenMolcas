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

!#define _DEBUGPRINT_
subroutine FIND_ACTIVE_BLOCKS(LUIN,LBLK,BLK_A,SEGMNT)
! Find the active (nonvanishing blocks) on LUIN
! Non vanishing block is flagged by a 1.0 (note : real)
! in BLK_A

use lucia_data, only: IDISK
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: LUIN, LBLK
real(kind=wp), intent(_OUT_) :: BLK_A(*), SEGMNT(*)
integer(kind=iwp) :: IAMPACK, IBLK, IDUMMY(1), IMZERO, KBLK, LBL(1), NBLK_A, NO_ZEROING
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: NBLK
#endif

IDISK(LUIN) = 0

IBLK = 0
NBLK_A = 0
! Loop over blocks
do
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
      BLK_A(IBLK) = One
    else
      BLK_A(IBLK) = Zero
    end if
  end if
  if ((LBL(1) < 0) .or. (LBLK > 0)) exit
end do

#ifdef _DEBUGPRINT_
NBLK = IBLK-1
write(u6,*) ' FIND_A.... Number of total and active Blocks',NBLK,NBLK_A
write(u6,*) ' Active blocks'
call WRTMAT(BLK_A,1,NBLK,1,NBLK)
#endif

end subroutine FIND_ACTIVE_BLOCKS
