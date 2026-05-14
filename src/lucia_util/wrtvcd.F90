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

subroutine WRTVCD(SEGMNT,LU,IREW,LBLK)
! PRINT VECTOR ON FILE LU
!
! LBLK DEFINES STRUCTURE OF FILES :

use lucia_data, only: IDISK
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: SEGMNT(*)
integer(kind=iwp), intent(in) :: LU, IREW, LBLK
integer(kind=iwp) :: IAMPACK, IBLK, IDUMMY(1), IMZERO, KBLK, LBL

if (IREW /= 0) then
  if (LBLK >= 0) then
    IDISK(LU) = 0
  else
    IDISK(LU) = 0
  end if
end if
! LOOP OVER BLOCKS

IBLK = 0
do
  if (LBLK > 0) then
    LBL = LBLK
  else if (LBLK == 0) then
    call IDAFILE(LU,2,IDUMMY,1,IDISK(LU))
    LBL = IDUMMY(1)
  else
    call IDAFILE(LU,2,IDUMMY,1,IDISK(LU))
    LBL = IDUMMY(1)
    call IDAFILE(LU,2,IDUMMY,1,IDISK(LU))
  end if
  IBLK = IBLK+1
  if (LBL >= 0) then
    if (LBLK >= 0) then
      KBLK = LBL
    else
      KBLK = -1
    end if
    call FRMDSC(SEGMNT,LBL,KBLK,LU,IMZERO,IAMPACK)
    if (LBL > 0) then
      write(u6,'(A,I3,A,I6)') ' Number of elements in segment ',IBLK,' IS ',LBL
      call WRTMAT(SEGMNT,1,LBL,1,LBL)
    end if
  end if

  if ((LBL < 0) .or. (LBLK > 0)) exit
end do

end subroutine WRTVCD
