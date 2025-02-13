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

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NSMST, NGRP, NELMNT, KGRP(NELMNT), ISMDFGP(NSMST,NGRP), NACTSYM(NGRP)
integer(kind=iwp), intent(inout) :: SCR(NELMNT)
integer(kind=iwp), intent(out) :: NONEW
integer(kind=iwp) :: I, IGAS, IPLACE, J, NTEST

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

if (NTEST /= 0) then
  write(u6,*) 'New ISMDFGP'
  write(u6,'(40I2)') (ISMDFGP(SCR(IGAS),KGRP(IGAS)),IGAS=1,NELMNT)
  write(u6,*) ' New SCR'
  call IWRTMA(SCR,1,NELMNT,1,NELMNT)
end if

end subroutine NXTDIST
