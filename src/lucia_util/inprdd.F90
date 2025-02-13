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

function INPRDD(VEC1,VEC2,LU1,LU2,IREW,LBLK)
! DISC VERSION OF INPROD
!
! LBLK DEFINES STRUCTURE OF FILE

use lucia_data, only: IDISK
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp) :: INPRDD
real(kind=wp), intent(_OUT_) :: VEC1(*), VEC2(*)
integer(kind=iwp), intent(in) :: LU1, LU2, IREW, LBLK
integer(kind=iwp) :: IAMPACK, IDUMMY(1), IMZERO, KBLK, NBL1
real(kind=wp) :: dDot_, X
logical(kind=iwp) :: DIFVEC

X = Zero
if (LU1 /= LU2) then
  DIFVEC = .true.
else
  DIFVEC = .false.
end if

if (IREW /= 0) then
  if (LBLK >= 0) then
    IDISK(LU1) = 0
    if (DIFVEC) IDISK(LU2) = 0
  else
    IDISK(LU1) = 0
    if (DIFVEC) IDISK(LU2) = 0
  end if
end if

! LOOP OVER BLOCKS OF VECTORS

do

  if (LBLK > 0) then
    NBL1 = LBLK
  else if (LBLK == 0) then
    call IDAFILE(LU1,2,IDUMMY,1,IDISK(LU1))
    NBL1 = IDUMMY(1)
    if (DIFVEC) call IDAFILE(LU2,2,IDUMMY,1,IDISK(LU2))
  else if (LBLK < 0) then
    call IDAFILE(LU1,2,IDUMMY,1,IDISK(LU1))
    NBL1 = IDUMMY(1)
    call IDAFILE(LU1,2,IDUMMY,1,IDISK(LU1))
    if (DIFVEC) then
      call IDAFILE(LU2,2,IDUMMY,1,IDISK(LU2))
      call IDAFILE(LU2,2,IDUMMY,1,IDISK(LU2))
    end if
  end if

  if (NBL1 >= 0) then
    if (LBLK >= 0) then
      KBLK = NBL1
    else
      KBLK = -1
    end if
    call FRMDSC(VEC1,NBL1,KBLK,LU1,IMZERO,IAMPACK)
    if (DIFVEC) then
      call FRMDSC(VEC2,NBL1,KBLK,LU2,IMZERO,IAMPACK)
      if (NBL1 > 0) X = X+dDot_(NBL1,VEC1,1,VEC2,1)
      !write(u6,*) ' vec1 and vec2 in INPRDD'
      !call WRTMAT(VEC1,1,NBL1,1,NBL1)
      !call WRTMAT(VEC2,1,NBL1,1,NBL1)
    else
      if (NBL1 > 0) X = X+dDot_(NBL1,VEC1,1,VEC1,1)
    end if
  end if
  if ((NBL1 < 0) .or. (LBLK > 0)) exit

end do

INPRDD = X

end function INPRDD
