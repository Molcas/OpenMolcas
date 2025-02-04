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

subroutine NXTDIST(NSMST,NGRP,NELMNT,KGRP,ISMDFGP,SCR,NACTSYM,NONEW)
! NONEW = 1 on return indicates that no additional numbers
! could be obtained.
!
!Giovanni Li Manni Feb 2012
!
!*************
! Input
!*************
! NSMST   = Number of Irreps
! NGRP    = Number of Groups
! NELMNT  = Number of GAS spaces
! KGRP    = K supergroup
! ISMDFGP = Symmetry distributions per group
!*************
! OUTPUT
!*************
! NONEW

use Definitions, only: u6

integer KGRP(NELMNT), ISMDFGP(NSMST,NGRP), NELMNT
integer SCR(NELMNT), NACTSYM(NGRP)

NTEST = 0
if (NTEST /= 0) then
  write(u6,*) 'NACTSYM :'
  call IWRTMA(NACTSYM,1,NGRP,1,NGRP)
  write(u6,*) ' ISMDFGP  Table:'
  do J=1,NSMST
    write(u6,'(40I2)') (ISMDFGP(J,I),I=1,NGRP)
  end do
  write(u6,*) 'Initial SCR'
  call IWRTMA(SCR,1,NELMNT,1,NELMNT)
end if

if (NELMNT == 0) then
  NONEW = 1
  goto 1001
end if

IPLACE = 0
1000 continue
IPLACE = IPLACE+1
if (SCR(IPLACE) < NACTSYM(KGRP(IPLACE))) then
  SCR(IPLACE) = SCR(IPLACE)+1
  NONEW = 0
  goto 1001
else if (IPLACE < NELMNT) then
  do JPLACE=1,IPLACE
    SCR(JPLACE) = 1
  end do
else if (IPLACE == NELMNT) then
  NONEW = 1
  goto 1001
end if
goto 1000
1001 continue

if (NTEST /= 0) then
  write(u6,*) 'New ISMDFGP'
  write(u6,'(40I2)') (ISMDFGP(SCR(IGAS),KGRP(IGAS)),IGAS=1,NELMNT)
end if

if (NTEST /= 0) then
  write(u6,*) ' New SCR'
  call IWRTMA(SCR,1,NELMNT,1,NELMNT)
end if

end subroutine NXTDIST
