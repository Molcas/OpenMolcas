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
! Copyright (C) 2003, Jeppe Olsen                                      *
!               2003, Jesper Wisborg Krogh                             *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine MXMNOC_OCCLS(MINEL,MAXEL,NORBTP,NORBFTP,NELFTP,MINOP)
! Construct accumulated MAX and MIN arrays for an occupation class
!
! MINOP (Smallest allowed number of open orbitals) added
! April2, 2003, JO (modified by JWK, April - June 2003)

use lucia_data, only: MXPNGAS
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(_OUT_) :: MINEL(*), MAXEL(*)
integer(kind=iwp), intent(in) :: NORBTP, NORBFTP(*), NELFTP(*), MINOP
integer(kind=iwp) :: IBORB, IGAS, IORB, IORB_START, MAX_DOUBLE, MAXOP_EXL, MAXOP_GAS(MXPNGAS), MAXOP_T, MINOP_GAS(MXPNGAS), &
                     NEL_INI, NELEC, NGAS
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: NORB

write(u6,*)
write(u6,*) ' ============'
write(u6,*) ' MXMNOC_OCCLS'
write(u6,*) ' ============'
write(u6,*)
write(u6,*) ' MINOP  = ',MINOP
write(u6,*) ' NORBTP = ',NORBTP
write(u6,*) ' NORBFTP :'
call IWRTMA(NORBFTP,1,NORBTP,1,NORBTP)
#endif
! Well
NGAS = NORBTP

! Largest number of unpaired electrons in each gas space

MAXOP_GAS(1:NGAS) = min(NELFTP(1:NGAS),2*NORBFTP(1:NGAS)-NELFTP(1:NGAS))

! Smallest number of electrons in each GAS space

! 1 : Just based on number of electrons in each space
MINOP_GAS(1:NGAS) = mod(NELFTP(1:NGAS),2)
! 2 : the total number of open orbitals should be MINOP, this puts
! also a constraint on the number of open orbitals

! The largest number of open orbitals, all spaces
MAXOP_T = sum(MAXOP_GAS(1:NGAS))
do IGAS=1,NGAS
  ! Max number of open orbitals in all spaces except IGAS
  MAXOP_EXL = MAXOP_T-MAXOP_GAS(IGAS)
  MINOP_GAS(IGAS) = max(MINOP_GAS(IGAS),MINOP-MAXOP_EXL)
  if (mod(NELFTP(IGAS)-MINOP_GAS(IGAS),2) == 1) MINOP_GAS(IGAS) = MINOP_GAS(IGAS)+1
end do
! We now have the min and max number of open shells per occls,
! Find the corresponding min and max number accumulated electrons,

! The Max occupation is obtained by occupying in max in the
! first orbitals
! The Min occupation is obtained by occopying max in the
! last orbitals.

NEL_INI = 0
IBORB = 1
do IGAS=1,NGAS
  NELEC = NELFTP(IGAS)
  ! PAM2009: This looks like a bug. MAX_DOUBLE can go negative.
  !MAX_DOUBLE = (NELEC-MINOP_GAS(IGAS))/2
  ! Replace with:
  MAX_DOUBLE = max(0,(NELEC-MINOP_GAS(IGAS))/2)

  ! If you are in a situation with no electrons to spare

  if (NELEC == 0) then
    do IORB=1,NORBFTP(IGAS)
      if (IORB+IBORB-1 == 1) then
        MINEL(IORB+IBORB-1) = 0
        MAXEL(IORB+IBORB-1) = 0
      else
        MINEL(IORB+IBORB-1) = MINEL(IORB+IBORB-2)
        MAXEL(IORB+IBORB-1) = MAXEL(IORB+IBORB-2)
      end if
    end do
  else
    ! The min number of electrons

    ! Doubly occupy the last MAX_DOUBLE orbitals
    ! Start Jesper !!!
    if ((NORBFTP(IGAS)-MAX_DOUBLE <= 0) .and. (MINOP_GAS(IGAS) > 0)) call Abend()
    ! End Jesper !!!
    IORB_START = max(1,NORBFTP(IGAS)-MAX_DOUBLE)
    do IORB=IORB_START,NORBFTP(IGAS)
      MINEL(IORB+IBORB-1) = NEL_INI+NELEC-2*(NORBFTP(IGAS)-IORB)
      !write(u6,*) ' 1 IORB+IBORB-1, MINEL() ',IORB+IBORB-1,MINEL(IORB+IBORB-1)
    end do
    ! Singly occupy
    do IORB=NORBFTP(IGAS)-MAX_DOUBLE-1,1,-1
      MINEL(IORB+IBORB-1) = max(NEL_INI,MINEL(IORB+IBORB-1+1)-1)
      !write(u6,*) ' 2 IORB+IBORB-1, MINEL() ',IORB+IBORB-1,MINEL(IORB+IBORB-1)
    end do

    ! The max number of electrons

    do IORB=1,MAX_DOUBLE
      MAXEL(IORB+IBORB-1) = NEL_INI+2*IORB
    end do
    do IORB=MAX_DOUBLE+1,NORBFTP(IGAS)
      if (IORB+IBORB-1 == 1) then
        MAXEL(IORB+IBORB-1) = 1
      else
        MAXEL(IORB+IBORB-1) = min(NEL_INI+NELEC,MAXEL(IORB+IBORB-2)+1)
      end if
    end do
  end if
  NEL_INI = NEL_INI+NELFTP(IGAS)
  IBORB = IBORB+NORBFTP(IGAS)
end do

#ifdef _DEBUGPRINT_
NORB = sum(NORBFTP(1:NORBTP))
write(u6,*) ' MINEL :'
call IWRTMA(MINEL,1,NORB,1,NORB)
write(u6,*) ' MAXEL :'
call IWRTMA(MAXEL,1,NORB,1,NORB)
#endif

end subroutine MXMNOC_OCCLS
