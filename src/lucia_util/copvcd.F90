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

subroutine COPVCD(LUIN,LUOUT,SEGMNT,IREW,LBLK)
! COPY VECTOR ON FILE LUIN TO FILE LUOUT
!
! LBLK DEFINES STRUCTURE OF FILE
!
! Structure of output file is inherited by output file,
! if input file is packed, so is output file
!
! Type of file LUOUT is inherited from LUIN

use lucia_data, only: IDISK
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: LUIN, LUOUT, IREW, LBLK
real(kind=wp), intent(_OUT_) :: SEGMNT(*)
integer(kind=iwp) :: IAMPACK, IDUMMY(1), IMZERO, KBLK, LBL(1), NO_ZEROING

if (IREW /= 0) then
  IDISK(LUIN) = 0
  IDISK(LUOUT) = 0
end if

! LOOP OVER BLOCKS

!write(u6,*) ' COPVCD LBLK : ',LBLK
do
  if (LBLK > 0) then
    LBL(1) = LBLK
  else if (LBLK == 0) then
    call IDAFILE(LUIN,2,LBL,1,IDISK(LUIN))
    call IDAFILE(LUOUT,1,LBL,1,IDISK(LUOUT))
    !write(u6,*) ' COPVCD LBL : ',LBL(1)
  else if (LBLK < 0) then
    call IDAFILE(LUIN,2,LBL,1,IDISK(LUIN))
    call IDAFILE(LUIN,2,IDUMMY,1,IDISK(LUIN))
    call IDAFILE(LUOUT,1,LBL,1,IDISK(LUOUT))
    IDUMMY(1) = -1
    call IDAFILE(LUOUT,1,IDUMMY,1,IDISK(LUOUT))
  end if
  if (LBL(1) >= 0) then
    if (LBLK >= 0) then
      KBLK = LBL(1)
    else
      KBLK = -1
    end if
    !write(u6,*) ' LBL and KBLK ',LBL(1),KBLK
    NO_ZEROING = 1
    call FRMDSC2(SEGMNT,LBL(1),KBLK,LUIN,IMZERO,IAMPACK,NO_ZEROING)
    !if (IAMPACK /= 0) write(u6,*) ' COPVCD, IAMPACK,FILE = ',IAMPACK,LUIN
    if (IMZERO == 0) then
      if (IAMPACK == 0) then
        call TODSC(SEGMNT,LBL(1),KBLK,LUOUT)
      else
        call TODSCP(SEGMNT,LBL(1),KBLK,LUOUT)
      end if
    else
      call ZERORC(LUOUT,IAMPACK)
    end if
  end if
  if ((LBL(1) < 0) .or. (LBLK > 0)) exit
end do

end subroutine COPVCD
