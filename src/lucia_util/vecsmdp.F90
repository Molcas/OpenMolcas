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

subroutine VECSMDP(VEC1,VEC2,FAC1,FAC2,LU1,LU2,LU3,IREW,LBLK)
! DISC VERSION OF VECSUM :
!
! ADD BLOCKED VECTORS ON FILES LU1 AND LU2
! AND STORE ON LU3
!
! Packed version, May 1996
!
! LBLK DEFINES STRUCTURE OF FILE

use lucia_data, only: IDISK
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: VEC1(*), VEC2(*)
real(kind=wp), intent(in) :: FAC1, FAC2
integer(kind=iwp), intent(in) :: LU1, LU2, LU3, IREW, LBLK
integer(kind=iwp) :: IAMPACK, IDUMMY(1), IMZERO1, IMZERO2, KBLK, NBL1, NBL2, NO_ZEROING

if (IREW /= 0) then
  IDISK(LU1) = 0
  IDISK(LU2) = 0
  IDISK(LU3) = 0
end if

! LOOP OVER BLOCKS OF VECTOR

do

  if (LBLK > 0) then
    NBL1 = LBLK
    NBL2 = LBLK
  else if (LBLK == 0) then
    call IDAFILE(LU1,2,IDUMMY,1,IDISK(LU1))
    NBL1 = IDUMMY(1)
    call IDAFILE(LU2,2,IDUMMY,1,IDISK(LU2))
    NBL2 = IDUMMY(1)
    IDUMMY(1) = NBL1
    call IDAFILE(LU3,1,IDUMMY,1,IDISK(LU3))
  else if (LBLK < 0) then
    call IDAFILE(LU1,2,IDUMMY,1,IDISK(LU1))
    NBL1 = IDUMMY(1)
    call IDAFILE(LU1,2,IDUMMY,1,IDISK(LU1))
    call IDAFILE(LU2,2,IDUMMY,1,IDISK(LU2))
    NBL2 = IDUMMY(1)
    call IDAFILE(LU2,2,IDUMMY,1,IDISK(LU2))
    IDUMMY(1) = NBL1
    call IDAFILE(LU3,1,IDUMMY,1,IDISK(LU3))
    IDUMMY(1) = -1
    call IDAFILE(LU3,1,IDUMMY,1,IDISK(LU3))
  end if
  if (NBL1 /= NBL2) then
    write(u6,'(A,2I5)') 'DIFFERENT BLOCKSIZES IN VECSMD ',NBL1,NBL2
    !stop ' INCOMPATIBLE BLOCKSIZES IN VECSMF'
    call SYSABENDMSG('lucia_util/vecsmf','Different block sizes','')
  end if

  if (NBL1 >= 0) then
    if (LBLK >= 0) then
      KBLK = NBL1
    else
      KBLK = -1
    end if
    NO_ZEROING = 1
    call FRMDSC2(VEC1,NBL1,KBLK,LU1,IMZERO1,IAMPACK,NO_ZEROING)
    call FRMDSC2(VEC2,NBL1,KBLK,LU2,IMZERO2,IAMPACK,NO_ZEROING)
    if (NBL1 > 0) then
      if ((IMZERO1 == 1) .and. (IMZERO2 == 1)) then
        ! Simple zero record
        call ZERORC(LU3,IAMPACK)
      else
        ! Nonvanishing record
        if (IMZERO1 == 1) then
          VEC1(1:NBL1) = FAC2*VEC2(1:NBL1)
        else if (IMZERO2 == 1) then
          VEC1(1:NBL1) = FAC1*VEC1(1:NBL1)
        else
          VEC1(1:NBL1) = FAC1*VEC1(1:NBL1)+FAC2*VEC2(1:NBL1)
        end if
        call TODSCP(VEC1,NBL1,KBLK,LU3)
      end if
    else if (NBL1 == 0) then
      call TODSCP(VEC1,NBL1,KBLK,LU3)
    end if
  end if

  if ((NBL1 < 0) .or. (LBLK > 0)) exit

end do

end subroutine VECSMDP
