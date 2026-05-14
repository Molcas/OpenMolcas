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
! Copyright (C) 1995, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine OCCLS(IWAY,NOCCLS,IOCCLS,NEL,NGAS,IGSMIN,IGSMAX,I_DO_BASSPC,IBASSPC,NOBPT)
! IWAY = 1 :
! obtain NOCCLS =
! Number of allowed ways of distributing the orbitals in the
! active spaces
!
! IWAY = 2 :
! OBTAIN NOCCLS and
! IOCCLS = allowed distributions of electrons
!
! Added Oct 98 : IBASSPC
! The basespace of
! a given class is the first space where this class occurs
!
! Jeppe Olsen, August 1995

use lucia_data, only: MXPNGAS
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: IWAY, NGAS, NEL, IGSMIN(NGAS), IGSMAX(NGAS), I_DO_BASSPC, NOBPT(NGAS)
integer(kind=iwp), intent(inout) :: NOCCLS, IBASSPC(*)
integer(kind=iwp), intent(_OUT_) :: IOCCLS(NGAS,*)
integer(kind=iwp) :: IBASSPC_FOR_CLS, IEL, IFIRST, IGAS, IM_TO_STUFFED, IOC(MXPNGAS), IOCA(MXPNGAS), ISKIP, KGAS, NEGA, NONEW
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: I

write(u6,*) ' OCCLS in action'
write(u6,*) ' =================='
write(u6,*) ' NGAS NEL ',NGAS,NEL
#endif

ISKIP = 1
NOCCLS = 0
! start with smallest allowed number
IOCA(1:NGAS) = IGSMIN(1:NGAS)
NONEW = 0
IFIRST = 1
! Loop over possible occupations
do
  if (IFIRST == 0) then
    ! Next accumulated occupation
    call NXTNUM3(IOCA,NGAS,IGSMIN,IGSMAX,NONEW)
  end if
  if (NONEW /= 0) exit
  ! ensure that IOCA corresponds to an accumulating occupation,
  ! i.e. a non-decreasing sequence
  if (ISKIP == 1) then
    KGAS = 0
    do IGAS=2,NGAS
      if (IOCA(IGAS-1) > IOCA(IGAS)) KGAS = IGAS
    end do
    if (KGAS /= 0) then
      IOCA(1:KGAS-1) = IGSMIN(1:KGAS-1)
      IOCA(KGAS) = IOCA(KGAS)+1
    end if
  end if
  !write(u6,*) ' Another accumulated occupation:'
  !call IWRTMA(IOCA,1,NGAS,1,NGAS)
  ! corresponding occupation of each active space
  NEGA = 0
  IM_TO_STUFFED = 0
  do IGAS=1,NGAS
    if (IGAS == 1) then
      IOC(IGAS) = IOCA(IGAS)
    else
      IOC(IGAS) = IOCA(IGAS)-IOCA(IGAS-1)
      if (IOC(IGAS) < 0) NEGA = 1
      if (IOC(IGAS) > 2*NOBPT(IGAS)) IM_TO_STUFFED = 1
    end if
  end do
  !write(u6,*) ' Another occupation:'
  !call IWRTMA(IOC,1,NGAS,1,NGAS)
  IFIRST = 0
  ! Correct number of electrons
  IEL = sum(IOC(1:NGAS))
  if ((IEL == NEL) .and. (NEGA == 0) .and. (IM_TO_STUFFED == 0)) then
    NOCCLS = NOCCLS+1
    if (IWAY == 2) then
#     ifdef _DEBUGPRINT_
      write(u6,*) ' Another allowed class :'
      call IWRTMA(IOC,1,NGAS,1,NGAS)
#     endif
      IOCCLS(:,NOCCLS) = IOC(1:NGAS)

      if (I_DO_BASSPC == 1) IBASSPC(NOCCLS) = IBASSPC_FOR_CLS(IOC)

    end if
  end if
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' Number of Allowed occupation classes ',NOCCLS
if (IWAY == 2) then
  write(u6,*) ' Occupation classes :'
  write(u6,*) ' ===================='
  write(u6,*)
  write(u6,*) ' Class    Occupation in GASpaces'
  write(u6,*) ' ================================'
  do I=1,NOCCLS
    write(u6,'(1X,I5,3X,16I3)') I,(IOCCLS(IGAS,I),IGAS=1,NGAS)
  end do
  !call IWRTMA(IOCCLS,NGAS,NOCCLS,NGAS,NOCCLS)
end if
#endif

!if (I_DO_BASSPC == 1) then
!  write(u6,*) ' Base CI spaces for the classes'
!  call IWRTMA(IBASSPC,1,NOCCLS,1,NOCCLS)
!end if

end subroutine OCCLS
