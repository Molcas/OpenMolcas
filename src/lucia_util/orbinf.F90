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

!#define _DEBUGPRINT_
subroutine ORBINF()
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

use lucia_data, only: IBSO, IOBPTS, IREOST, IREOTS, ISMFSO, ISMFTO, ITOOBS, ITOOBS, MXPIRR, MXPNGAS, MXPOBS, MXTSOB, NACOB, &
                      NACOBS, NDEOB, NDEOB, NGAS, NGSOBT, NGSSH, NINOB, NINOB, NINOBS, NIRREP, NOBPT, NOBPTS, NOCOB, NSMOB, NTOOB, &
                      NTOOBS
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: IGAS, IGSDEL, IGSINA, IOBSM, IOBTP, IOSPIR(MXPOBS,MXPIRR), ISMOB, LTOB, NGSOB(MXPOBS,MXPNGAS), NOSPIR(MXPIRR)
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: I
#endif

!***********************************************
!                                              *
! Part 1 : From shell format to orbital format *
!                                              *
!***********************************************
call OSPIR(NOSPIR,IOSPIR,MXPIRR,MXPOBS)

! 2 : Shell information to orbital information for each group of orbital

! ===============
!     GAS case
! ===============

do IGAS=1,NGAS
  ! Shell => orbitals for each GAS space
  call SHTOOB(NGSSH(:,IGAS),NIRREP,MXPOBS,NSMOB,NOSPIR,IOSPIR,NGSOB(:,IGAS),NGSOBT(IGAS))
end do

! ======================================================
! Number of inactive, active, occupied, deleted orbitals
! ======================================================

! current inactive and deleted orbitals are not identified so
IGSINA = 0
IGSDEL = 0

NTOOBS(1:NSMOB) = 0
!NOCOBS(1:NSMOB) = 0
NACOBS(1:NSMOB) = 0

NTOOB = 0
NACOB = 0
NOCOB = 0
do IGAS=1,NGAS
  ! Inactive orbitals
  if (IGAS == IGSINA) then
    NINOBS(1:NSMOB) = NGSOB(1:NSMOB,IGAS)
    NINOB = NGSOBT(IGAS)
  end if
  ! Deleted orbitals
  if (IGAS == IGSDEL) then
    !NDEOBS(1:NSMOB) = NGSOB(1:NSMOB,IGAS)
    NDEOB = NGSOBT(IGAS)
  end if
  ! Add to total number of orbitals
  NTOOBS(1:NSMOB) = NTOOBS(1:NSMOB)+NGSOB(1:NSMOB,IGAS)
  NTOOB = NTOOB+NGSOBT(IGAS)
  ! Add to occupied orbitals
  if (IGAS /= IGSDEL) then
    !NOCOBS(1:NSMOB) = NOCOBS(1:NSMOB)+NGSOB(1:NSMOB,IGAS)
    NOCOB = NOCOB+NGSOBT(IGAS)
  end if
  ! Add to active orbitals
  if ((IGAS /= IGSINA) .and. (IGAS /= IGSDEL)) then
    NACOBS(1:NSMOB) = NACOBS(1:NSMOB)+NGSOB(1:NSMOB,IGAS)
    NACOB = NACOB+NGSOBT(IGAS)
  end if
end do
#ifdef _DEBUGPRINT_
! =================
! Well, report back
! =================
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
#endif
! Offsets for orbitals of given symmetry
ITOOBS(1) = 1
do ISMOB=2,NSMOB
  ITOOBS(ISMOB) = ITOOBS(ISMOB-1)+NTOOBS(ISMOB-1)
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' Offsets for orbital of given symmetry'
call IWRTMA(ITOOBS,1,NSMOB,1,NSMOB)
#endif

!*******************************************
!                                          *
! Part 2 : Reordering arrays for orbitals  *
!                                          *
!*******************************************
call ORBORD_GAS(NSMOB,MXPOBS,MXPNGAS,NGAS,NGSOB,NGSOBT,NTOOBS,NTOOB,IREOST,IREOTS,ISMFTO,IBSO,NOBPTS,IOBPTS,ISMFSO,NOBPT)

! Largest number of orbitals of given sym and type
MXTSOB = 0
do IOBTP=1,NGAS
  LTOB = 0
  do IOBSM=1,NSMOB
    MXTSOB = max(MXTSOB,NOBPTS(IOBTP,IOBSM))
    LTOB = LTOB+NOBPTS(IOBTP,IOBSM)
  end do
end do
#ifdef _DEBUGPRINT_
write(u6,*) ' MXTSOB from ORBINF = ',MXTSOB
#endif

end subroutine ORBINF
