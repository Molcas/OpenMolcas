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

!#define _DEBUGPRINT_
subroutine NXTNUM3(INUM,NELMNT,MNVAL,MXVAL,NONEW)
! A set of numbers INUM(I),I=1,NELMNT is given.
! Find next compund number.
! Digit I must be in the range MNVAL(I),MXVAL(I).
!
! NONEW = 1 on return indicates that no additional numbers
! could be obtained.
!
! Jeppe Olsen Oct 1994

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NELMNT, MNVAL(NELMNT), MXVAL(NELMNT)
integer(kind=iwp), intent(inout) :: INUM(NELMNT)
integer(kind=iwp), intent(out) :: NONEW
integer(kind=iwp) :: IPLACE

#ifdef _DEBUGPRINT_
write(u6,*) ' Initial number in NXTNUM'
call IWRTMA(INUM,1,NELMNT,1,NELMNT)
#endif

if (NELMNT == 0) then
  NONEW = 1
else

  IPLACE = 0
  do
    IPLACE = IPLACE+1
    if (INUM(IPLACE) < MXVAL(IPLACE)) then
      INUM(IPLACE) = INUM(IPLACE)+1
      NONEW = 0
      exit
    else if (IPLACE < NELMNT) then
      INUM(1:IPLACE) = MNVAL(1:IPLACE)
    else if (IPLACE == NELMNT) then
      NONEW = 1
      exit
    end if
  end do

end if

#ifdef _DEBUGPRINT_
write(u6,*) ' New number'
call IWRTMA(INUM,1,NELMNT,1,NELMNT)
#endif

end subroutine NXTNUM3
