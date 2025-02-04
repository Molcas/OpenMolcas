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
! Copyright (C) 1991, Jeppe Olsen                                      *
!***********************************************************************

subroutine ORBINF(IPRNT)
! Obtain information about orbitals from shell information
!
! =====
! input
! =====
! Shell and symmetry information in /LUCINP/
!
! ======
! Output
! ======
! Orbital information in /ORBINP/
!
! Jeppe Olsen, Winter of 1991

use lucia_data, only: NGAS, IGSINA, IGSDEL, NGSOB, NGSOBT, NGSSH
use lucia_data, only: NSMOB, NIRREP, PNTGRP
use lucia_data, only: NTOOB, NACOB, NOCOB, NINOB, NDEOB, MXTSOB, MXTOB, ITOOBS, IBSO, IOBPTS, IOSPIR, IREOST, IREOTS, ISMFSO, &
                      ISMFTO, ITOOBS, ITPFSO, ITPFTO, NACOBS, NDEOB, NDEOBS, NINOB, NINOBS, NOBPT, NOBPTS, NOCOBS, NOSPIR, NTOOBS
use lucia_data, only: MXPIRR, MXPNGAS, MXPOBS
use Definitions, only: u6

implicit none
integer IPRNT
integer NTEST, IGAS, I, ISMOB, IOBTP, LTOB, IOBSM

NTEST = 0
NTEST = max(NTEST,IPRNT)
!***********************************************
!                                              *
! Part 1 : From shell format to orbital format *
!                                              *
!***********************************************
call OSPIR(NOSPIR,IOSPIR,PNTGRP,NIRREP,MXPIRR,MXPOBS,IPRNT)

! 2 : Shell information to orbital information for each group of orbital

! ===============
!     GAS case
! ===============

do IGAS=1,NGAS
  ! Shell => orbitals for each GAS space
  call SHTOOB(NGSSH(1,IGAS),NIRREP,MXPOBS,NSMOB,NOSPIR,IOSPIR,NGSOB(1,IGAS),NGSOBT(IGAS))
end do

! ======================================================
! Number of inactive, active, occupied, deleted orbitals
! ======================================================

! current inactive and deleted orbitals are not identified so
IGSINA = 0
IGSDEL = 0

call ISETVC(NTOOBS,0,NSMOB)
call ISETVC(NOCOBS,0,NSMOB)
call ISETVC(NACOBS,0,NSMOB)

NTOOB = 0
NACOB = 0
NOCOB = 0
do IGAS=1,NGAS
  ! Inactive orbitals
  if (IGAS == IGSINA) then
    call ICOPVE(NGSOB(1,IGAS),NINOBS,NSMOB)
    NINOB = NGSOBT(IGAS)
  end if
  ! Deleted orbitals
  if (IGAS == IGSDEL) then
    call ICOPVE(NGSOB(1,IGAS),NDEOBS,NSMOB)
    NDEOB = NGSOBT(IGAS)
  end if
  ! Add to total number of orbitals
  call IVCSUM(NTOOBS,NTOOBS,NGSOB(1,IGAS),1,1,NSMOB)
  NTOOB = NTOOB+NGSOBT(IGAS)
  ! Add to occupied orbitals
  if (IGAS /= IGSDEL) then
    call IVCSUM(NOCOBS,NOCOBS,NGSOB(1,IGAS),1,1,NSMOB)
    NOCOB = NOCOB+NGSOBT(IGAS)
  end if
  ! Add to active orbitals
  if ((IGAS /= IGSINA) .and. (IGAS /= IGSDEL)) then
    call IVCSUM(NACOBS,NACOBS,NGSOB(1,IGAS),1,1,NSMOB)
    NACOB = NACOB+NGSOBT(IGAS)
  end if
end do
! =================
! Well, report back
! =================
if (NTEST > 0) then
  write(u6,*)
  write(u6,*) ' Number of orbitals per symmetry :'
  write(u6,*) ' ================================='
  write(u6,*)
  write(u6,'(1X,A,10I4,A)') '            Symmetry  ',(I,I=1,NSMOB)
  write(u6,'(1X,A,2X,10A)') '           ========== ',('====',I=1,NSMOB)
  do IGAS=1,NGAS
    write(u6,'(1X,A,I3,7X,A,10I4,8X,I3)') '   GAS',IGAS,'      ',(NGSOB(I,IGAS),I=1,NSMOB),NGSOBT(IGAS)
  end do

  write(u6,*) ' Total number of orbitals ',NTOOB
  write(u6,*) ' Total number of occupied orbitals ',NOCOB
end if
! Offsets for orbitals of given symmetry
ITOOBS(1) = 1
do ISMOB=2,NSMOB
  ITOOBS(ISMOB) = ITOOBS(ISMOB-1)+NTOOBS(ISMOB-1)
end do

if (NTEST > 0) then
  write(u6,*) ' Offsets for orbital of given symmetry'
  call IWRTMA(ITOOBS,1,NSMOB,1,NSMOB)
end if

!*******************************************
!                                          *
! Part 2 : Reordering arrays for orbitals  *
!                                          *
!*******************************************
call ORBORD_GAS(NSMOB,MXPOBS,MXPNGAS,NGAS,NGSOB,NGSOBT,NOCOBS,NTOOBS,NTOOB,IREOST,IREOTS,ISMFTO,ITPFSO,IBSO,NOBPTS,IOBPTS,ISMFSO, &
                ITPFTO,NOBPT,IPRNT)

! Largest number of orbitals of given sym and type
MXTSOB = 0
MXTOB = 0
do IOBTP=1,NGAS
  LTOB = 0
  do IOBSM=1,NSMOB
    MXTSOB = max(MXTSOB,NOBPTS(IOBTP,IOBSM))
    LTOB = LTOB+NOBPTS(IOBTP,IOBSM)
  end do
  MXTOB = max(LTOB,MXTOB)
end do
if (NTEST > 0) write(u6,*) ' MXTSOB,MXTOB from ORBINF = ',MXTSOB,MXTOB

end subroutine ORBINF
