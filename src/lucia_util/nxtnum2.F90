!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1994, Jeppe Olsen                                      *
!***********************************************************************

subroutine NXTNUM2(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
! An set of numbers INUM(I),I=1,NELMNT is
! given. Find next compund number.
! Digit I must be in the range MINVAL,MAXVAL(I).
!
! NONEW = 1 on return indicates that no additional numbers
! could be obtained.
!
! Jeppe Olsen Oct 1994

use Definitions, only: u6

! Input
dimension maxval(*)
! Input and output
dimension INUM(*)

NTEST = 0
if (NTEST /= 0) then
  write(u6,*) ' Initial number in NXTNUM'
  call IWRTMA(INUM,1,NELMNT,1,NELMNT)
end if

if (NELMNT == 0) then
  NONEW = 1
  goto 1001
end if

IPLACE = 0
1000 continue
IPLACE = IPLACE+1
if (INUM(IPLACE) < maxval(IPLACE)) then
  INUM(IPLACE) = INUM(IPLACE)+1
  NONEW = 0
  goto 1001
else if (IPLACE < NELMNT) then
  do JPLACE=1,IPLACE
    !08-08-01 INUM(JPLACE) = 1
    INUM(JPLACE) = MINVAL
  end do
else if (IPLACE == NELMNT) then
  NONEW = 1
  goto 1001
end if
goto 1000
1001 continue

if (NTEST /= 0) then
  write(u6,*) ' New number'
  call IWRTMA(INUM,1,NELMNT,1,NELMNT)
end if

end subroutine NXTNUM2
