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

subroutine ORBORD(NSM,NR4TP,NRSOBS,NTOOBS,IREOTS,ISFTO,IBSO,NTSOB,IBTSOB,ITSOB,NOBPTS,NNOBPT,NOBPT)
! Obtain Reordering arrays for orbitals
! (See note below for assumed ordering)
!
! =====
! Input
! =====
!  NSM    : Number of orbital symmetries
!  NR4TP  : Number of RAS4 types
!  NRSOBS : Number of orbitals per symmetry in RAS1,RAS2,RAS3
!  NTOOBS : Number of orbitals per symmetry,all types
!
! ======
! Output
! ======
!  IREOTS : Reordering array type     => symmetry
!  ISFTO  : Symmetry array for type ordered orbitals
!  IBSO   : First orbital of given symmetry (symmetry ordered)
!  NTSOB  : Number of active orbitals of give RAS type and SYM
!  IBTSOB : Off set for active orbitals of given RAS type and SYM
!  ITSOB  : Orbitals of given RAS type and sym
!
!  NOBPTS : Number of orbitals per subtype and symmetry
!
! Jeppe Olsen, Winter 1991

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NSM, NR4TP, NRSOBS(NSM,3), NTOOBS(*), NNOBPT
integer(kind=iwp), intent(_OUT_) :: IREOTS(*), ISFTO(*), IBSO(*), NTSOB(3,*), IBTSOB(3,*), ITSOB(*), NOBPTS(NNOBPT,*)
integer(kind=iwp), intent(out) :: NOBPT(NNOBPT)
integer(kind=iwp) :: i, I123, IAC, IACS, IBSM, IIAC, IOFF, IORB, IOTYPE, IRS, ISM, ISMOB, ITYPE, LORB, NPREVS
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: NACOB, NTOOB
#endif

! =========================
! Note on order of orbitals
! =========================
!
! The orbitals are supposed to be imported ordered symmetry-type
! ordered as
! Loop over symmetries of orbitals
!   Inactive  of this symmetry
!   Core      of this symmetry
!   RAS1      of this symmetry
!   RAS2      of this symmetry
!   RAS3      of this symmetry
!   Secondary of this symmetry
!   Deleted   of this symmetry
! End of loop over symmetries
!
! Internally the orbitals are reordered to type symmetry order
! where the outer loop os over types and the inner loop is
! over symmetries.The types are arranged as
!  Ras1
!  Ras2
!  Ras3
!  Core
!  Secondary
!  Inactive
!  Deleted orbitals

! Active orbitals

IAC = 0
!IBSM = 0   ! dummy initialize
!NPREVS = 0 ! dummy initialize
!IORB = 0   ! dummy initialize
do IRS=1,3
  do ISM=1,NSM
    if (ISM == 1) then
      IBSM = 1
    else
      IBSM = IBSM+NTOOBS(ISM-1)
    end if
    NPREVS = sum(NRSOBS(ISM,1:IRS-1)) ! +NINOBS(ISM)+NR0OBS(ISM)
    IORB = NRSOBS(ISM,IRS)
    !if (IRS == 1) then
    !  NPREVS = 0 ! NINOBS(ISM)+NR0OBS(ISM)
    !  IORB = NRSOBS(ISM,1)
    !else if (IRS == 2) then
    !  NPREVS = NRSOBS(ISM,1) ! NINOBS(ISM)+NR0OBS(ISM)+NRSOBS(ISM,1)
    !  IORB = NRSOBS(ISM,2)
    !else if (IRS == 3) then
    !  NPREVS = NRSOBS(ISM,1)+NRSOBS(ISM,2) ! NINOBS(ISM)+NR0OBS(ISM)+NRSOBS(ISM,1)+NRSOBS(ISM,2)
    !  IORB = NRSOBS(ISM,3)
    !end if
    do IIAC=1,IORB
      ! Type ordered index
      IAC = IAC+1
      ! Symmetry ordered index
      IACS = IBSM+NPREVS-1+IIAC
      ISFTO(IAC) = ISM
      !ITFSO(IACS) = IRS
      IREOTS(IAC) = IACS
    end do
  end do
end do
!NACOB = IAC
!write(u6,*) ' IAC ',IAC

! RAS 0 orbitals

!IR0 = NACOB
!do ISM=1,NSM
!  if (ISM == 1) then
!    IBSM = 1
!  else
!    IBSM = IBSM+NTOOBS(ISM-1)
!  end if
!  NPREVS = NINOBS(ISM)
!  do IIR0=1,NR0OBS(ISM)
!    ! Type ordered index
!    IR0 = IR0+1
!    ! Symmetry ordered index
!    IR0S = IBSM+NPREVS-1+IIR0
!    ISFTO(IR0) = ISM
!    !ITFSO(IR0S) = 4
!    IREOTS(IR0) = IR0S
!  end do
!end do
!NR0OB = IR0-NACOB
!write(u6,*) ' IR0 ',IR0

! RAS 4 orbitals

!IR4 = NACOB+NR0OB
!do ISM=1,NSM
!  if (ISM == 1) then
!    IBSM = 1
!  else
!    IBSM = IBSM+NTOOBS(ISM-1)
!  end if
!  NPREVS = NINOBS(ISM)+NR0OBS(ISM)+NACOBS(ISM)
!  do ITP=1,NR4TP
!    do IIR4=1,NR4OBS(ISM,ITP)
!      ! Type ordered index
!      IR4 = IR4+1
!      ! Symmetry ordered index
!      IR4S = IBSM+NPREVS-1+IIR4
!      ISFTO(IR4) = ISM
!      !ITFSO(IR4S) = 5
!      IREOTS(IR4) = IR4S
!    end do
!  end do
!end do
!NR4OB = IR4-NACOB-NR0OB
!write(u6,*) ' IR4 ',IR4

! Inactive orbitals

!IIN = NACOB+NR0OB+NR4OB
!do ISM=1,NSM
!  if (ISM == 1) then
!    IBSM = 1
!  else
!    IBSM = IBSM+NTOOBS(ISM-1)
!  end if
!  NPREVS = 0
!  do IIIN=1,NINOBS(ISM)
!    ! Type ordered index
!    IIN = IIN+1
!    ! Symmetry ordered index
!    IINS = IBSM+NPREVS-1+IIIN
!    ISFTO(IIN) = ISM
!    !ITFSO(IINS) = 6
!    IREOTS(IIN) = IINS
!  end do
!end do
!NINOB = IIN-NACOB-NR0OB-NR4OB
!write(u6,*) ' IIN ',IIN

! Deleted orbitals

!IDE = NACOB+NR0OB+NR4OB+NINOB
!do ISM=1,NSM
!  if (ISM == 1) then
!    IBSM = 1
!  else
!    IBSM = IBSM+NTOOBS(ISM-1)
!  end if
!  IR4 = 0
!  do ITP=1,NR4TP
!    IR4 = IR4+NR4OBS(ISM,ITP)
!  end do
!
!  NPREVS = NINOBS(ISM)+NR0OBS(ISM)+NACOBS(ISM)+IR4
!  do IIDE=1,NDEOBS(ISM)
!    ! Type ordered index
!    IDE = IDE+1
!    ! Symmetry ordered index
!    IDES = IBSM+NPREVS-1+IIDE
!    ISFTO(IDE) = ISM
!    !ITFSO(IDES) = 7
!    IREOTS(IDE) = IDES
!  end do
!end do
!write(u6,*) ' IDE ',IDE

IBSO(1) = 1
do ISM=2,NSM
  IBSO(ISM) = IBSO(ISM-1)+NTOOBS(ISM-1)
end do

! ==================
! NTSOB,IBTSOB,ITSOB
! ==================

IOFF = 1
do I123=1,3
  do ISM=1,NSM
    NTSOB(I123,ISM) = NRSOBS(ISM,I123)
    IBTSOB(I123,ISM) = IOFF
    ITSOB(IOFF:IOFF+NRSOBS(ISM,I123)-1) = [(i,i=IOFF,IOFF+NRSOBS(ISM,I123)-1)]
    IOFF = IOFF+NRSOBS(ISM,I123)
  end do
end do

! ============
! NOBPTS NOBPT
! ============

! Loop over types in input order
NOBPT(1:NR4TP+6) = 0
LORB = 0   ! dummy initialize
IOTYPE = 0 ! dummy initialize
do ISMOB=1,NSM
  do ITYPE=1,NR4TP+6
    if (ITYPE == 1) then
      ! Inactive (frozen in normal notation)
      LORB = 0 ! NINOBS(ISMOB)
      IOTYPE = 5+NR4TP
    else if (ITYPE == 2) then
      ! RAS0 (inactive in normal notation)
      LORB = 0 ! NR0OBS(ISMOB)
      IOTYPE = 4
    else if (ITYPE == 3) then
      ! RAS1
      LORB = NRSOBS(ISMOB,1)
      IOTYPE = 1
    else if (ITYPE == 4) then
      ! RAS2
      LORB = NRSOBS(ISMOB,2)
      IOTYPE = 2
    else if (ITYPE == 5) then
      ! RAS3
      LORB = NRSOBS(ISMOB,3)
      IOTYPE = 3
    else if ((ITYPE >= 6) .and. (ITYPE <= NR4TP+6-1)) then
      ! RAS4
      LORB = 0 ! NR4OBS(ISMOB,ITYPE-5)
      IOTYPE = ITYPE-1
    else if (ITYPE == NR4TP+6) then
      ! deleted orbitals
      LORB = 0 ! NDEOBS(ISMOB)
      IOTYPE = ITYPE
    end if
    NOBPTS(IOTYPE,ISMOB) = LORB
    NOBPT(IOTYPE) = NOBPT(IOTYPE)+LORB
  end do
end do

#ifdef _DEBUGPRINT_
NACOB = IAC
NTOOB = NACOB
!NTOOB = IDE
write(u6,*) ' =================='
write(u6,*) ' Output from ORBORD'
write(u6,*) ' =================='
write(u6,*) ' Symmetry of orbitals, type ordered'
call IWRTMA(ISFTO,1,NTOOB,1,NTOOB)
write(u6,*) ' Type => symmetry reordering array'
call IWRTMA(IREOTS,1,NTOOB,1,NTOOB)
write(u6,*) ' IBSO array'
call IWRTMA(IBSO,1,NSM,1,NSM)

write(u6,*) ' NTSOB array :'
call IWRTMA(NTSOB,3,NSM,3,NSM)
write(u6,*) ' IBTSOB array'
call IWRTMA(IBTSOB,3,NSM,3,NSM)
write(u6,*) ' ITSOB'
call IWRTMA(ITSOB,1,NACOB,1,NACOB)

write(u6,*) ' NOBPTS'
call IWRTMA(NOBPTS,NR4TP+6,NSM,NNOBPT,NSM)
write(u6,*) ' NOBPT'
call IWRTMA(NOBPT,NR4TP+6,1,NNOBPT,1)

write(u6,*) ' ISFTO array :'
call IWRTMA(ISFTO,1,NTOOB,1,NTOOB)
!write(u6,*) ' ITFSO array :'
!call IWRTMA(ITFSO,1,NTOOB,1,NTOOB)
#endif

end subroutine ORBORD
