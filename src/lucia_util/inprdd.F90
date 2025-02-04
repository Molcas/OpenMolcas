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

real*8 function INPRDD(VEC1,VEC2,LU1,LU2,IREW,LBLK)
! DISC VERSION OF INPROD
!
! LBLK DEFINES STRUCTURE OF FILE

use lucia_data, only: IDISK
use Constants, only: Zero

implicit none
real*8 VEC1(*), VEC2(*)
integer LU1, LU2, IREW, LBLK
real*8 INPROD, X
logical DIFVEC
integer IDUMMY(1), NBL1, KBLK, IAMPACK, IMZERO

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

1000 continue

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
    if (NBL1 > 0) X = X+INPROD(VEC1,VEC2,NBL1)
    !write(u6,*) ' vec1 and vec2 in INPRDD'
    !call WRTMAT(VEC1,1,NBL1,1,NBL1)
    !call WRTMAT(VEC2,1,NBL1,1,NBL1)
  else
    if (NBL1 > 0) X = X+INPROD(VEC1,VEC1,NBL1)
  end if
end if
if ((NBL1 >= 0) .and. (LBLK <= 0)) goto 1000

INPRDD = X

end function INPRDD
