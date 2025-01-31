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
! Copyright (C) 1989, Jeppe Olsen                                      *
!***********************************************************************

subroutine NXTORD(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
! An ordered set of numbers INUM(I),I=1,NELMNT is
! given in strictly ascending order. Values of INUM(*) is
! restricted to the interval MINVAL,MAXVAL .
!
! Find next higher number.
!
! NONEW = 1 on return indicates that no additional numbers
! could be obtained.
!
! Jeppe Olsen May 1989

dimension INUM(*)

NTEST = 0
if (NTEST /= 0) then
  write(6,*) ' Initial number in NXTORD'
  call IWRTMA(INUM,1,NELMNT,1,NELMNT)
end if

IPLACE = 0
1000 continue
IPLACE = IPLACE+1
if ((IPLACE < NELMNT) .and. (INUM(IPLACE)+1 < INUM(IPLACE+1)) .or. (IPLACE == NELMNT) .and. (INUM(IPLACE)+1 <= MAXVAL)) then
  INUM(IPLACE) = INUM(IPLACE)+1
  NONEW = 0
  goto 1001
else if (IPLACE < NELMNT) then
  if (IPLACE == 1) then
    INUM(IPLACE) = MINVAL
  else
    INUM(IPLACE) = INUM(IPLACE-1)+1
  end if
else if (IPLACE == NELMNT) then
  NONEW = 1
  goto 1001
end if
goto 1000
1001 continue

if (NTEST /= 0) then
  write(6,*) ' New number'
  call IWRTMA(INUM,1,NELMNT,1,NELMNT)
end if

end subroutine NXTORD
