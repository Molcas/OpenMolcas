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
! Copyright (C) 2012, Giovanni Li Manni                                *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine NXTDIST(NGRP,NELMNT,KGRP,SCR,NACTSYM,NONEW)
! NONEW = 1 on return indicates that no additional numbers
! could be obtained.
!
!Giovanni Li Manni Feb 2012
!
!*************
! Input
!*************
! NGRP    = Number of Groups
! NELMNT  = Number of GAS spaces
! KGRP    = K supergroup
!*************
! OUTPUT
!*************
! NONEW

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NGRP, NELMNT, KGRP(NELMNT), NACTSYM(NGRP)
integer(kind=iwp), intent(inout) :: SCR(NELMNT)
integer(kind=iwp), intent(out) :: NONEW
integer(kind=iwp) :: IPLACE

#ifdef _DEBUGPRINT_
write(u6,*) 'NACTSYM :'
call IWRTMA(NACTSYM,1,NGRP,1,NGRP)
write(u6,*) 'Initial SCR'
call IWRTMA(SCR,1,NELMNT,1,NELMNT)
#endif

if (NELMNT == 0) then
  NONEW = 1
else

  IPLACE = 0
  do
    IPLACE = IPLACE+1
    if (SCR(IPLACE) < NACTSYM(KGRP(IPLACE))) then
      SCR(IPLACE) = SCR(IPLACE)+1
      NONEW = 0
      exit
    else if (IPLACE < NELMNT) then
      SCR(1:IPLACE) = 1
    else if (IPLACE == NELMNT) then
      NONEW = 1
      exit
    end if
  end do

end if

#ifdef _DEBUGPRINT_
write(u6,*) ' New SCR'
call IWRTMA(SCR,1,NELMNT,1,NELMNT)
#endif

end subroutine NXTDIST
