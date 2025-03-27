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

subroutine ORBORD(NSMOB,MXPOBS,NR4TP,NDEOBS,NINOBS,NR0OBS,NACOBS,NRSOBS,NR4OBS,NOCOBS,NTOOBS,IREOST,IREOTS,ISFTO,ITFSO,IBSO,NTSOB, &
                  IBTSOB,ITSOB,NOBPTS,IOBPTS,MXPR4T,ISMFSO,ITPFTO,NOBPT)
! Obtain Reordering arrays for orbitals
! (See note below for assumed ordering)
!
! =====
! Input
! =====
!  NSMOB  : Number of orbital symmetries
!  MXPOBS : Max number of orbital symmetries
!  NR4TP  : Number of RAS4 types
!  NDEOBS : Number of deleted orbitals per symmetry
!  NINOBS : Number of inactive  orbitals per symmetry
!  NR0OBS : Number of RAS 0 (core) orbitals per symmetry
!  NACOBS : Number of Active orbitals per symmetry
!  NRSOBS : Number of orbitals per symmetry in RAS1,RAS2,RAS3
!  NR4OBS : Number of RAS 4 orbitals per symmetry and type
!  NOCOBS : Number of occupied orbitals per symmetry
!  NTOOBS : Number of orbitals per symmetry,all types
!
! ======
! Output
! ======
!  IREOST : Reordering array symmetry => type
!  IREOTS : Reordering array type     => symmetry
!  ISFTO  : Symmetry array for type ordered orbitals
!  ITFSO  : Type array for symmetry ordered orbitals (not activated)
!  IBSO   : First orbital of given symmetry (symmetry ordered)
!  NTSOB  : Number of active orbitals of give RAS type and SYM
!  IBTSOB : Off set for active orbitals of given RAS type and SYM
!  ITSOB  : Orbitals of given RAS type and sym
!
!  NOBPTS : Number of orbitals per subtype and symmetry
!  IOBPTS : Off sets for orbitals of given subtype and symmetry
!           ordered according to input integrals
!
! ISMFSO  : Symmetry of orbitals, symmetry ordereing
! ITPFTO  : Type of orbital, type ordering
!
! Jeppe Olsen, Winter 1991

#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit real*8(A-H,O-Z)
! Input
dimension NDEOBS(*), NINOBS(*), NR0OBS(*), NACOBS(*), NRSOBS(MXPOBS,3), NR4OBS(MXPOBS,*), NOCOBS(*), NTOOBS(*)
! Output
dimension IREOST(*), IREOTS(*), ISFTO(*), ITFSO(*), IBSO(*)
dimension ISMFSO(*), ITPFTO(*)
dimension NTSOB(3,*), IBTSOB(3,*), ITSOB(*)
dimension NOBPTS(6+MXPR4T,*), IOBPTS(6+MXPR4T,*)
dimension NOBPT(6+MXPR4T)

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
IBSM = 0   ! dummy intitialize
NPREVS = 0 ! dummy initialize
IORB = 0   ! dummy initialize
do IRS=1,3
  do ISM=1,NSMOB
    if (ISM == 1) then
      IBSM = 1
    else
      IBSM = IBSM+NTOOBS(ISM-1)
    end if
    if (IRS == 1) then
      NPREVS = NINOBS(ISM)+NR0OBS(ISM)
      IORB = NRSOBS(ISM,1)
    else if (IRS == 2) then
      NPREVS = NINOBS(ISM)+NR0OBS(ISM)+NRSOBS(ISM,1)
      IORB = NRSOBS(ISM,2)
    else if (IRS == 3) then
      NPREVS = NINOBS(ISM)+NR0OBS(ISM)+NRSOBS(ISM,1)+NRSOBS(ISM,2)
      IORB = NRSOBS(ISM,3)
    end if
    do IIAC=1,IORB
      ! Type ordered index
      IAC = IAC+1
      ! Symmetry ordered index
      IACS = IBSM+NPREVS-1+IIAC
      ISFTO(IAC) = ISM
      !ISMFSO(IACS) = ISM
      ITPFTO(IAC) = IRS
      !ITFSO(IACS) = IRS
      IREOST(IACS) = IAC
      IREOTS(IAC) = IACS

    end do
  end do
end do
NACOB = IAC
!write(u6,*) ' IAC ',IAC

! RAS 0 orbitals

IR0 = NACOB
do ISM=1,NSMOB
  if (ISM == 1) then
    IBSM = 1
  else
    IBSM = IBSM+NTOOBS(ISM-1)
  end if
  NPREVS = NINOBS(ISM)
  do IIR0=1,NR0OBS(ISM)
    ! Type ordered index
    IR0 = IR0+1
    ! Symmetry ordered index
    IR0S = IBSM+NPREVS-1+IIR0
    ISFTO(IR0) = ISM
    !ITFSO(IR0S) = 4
    IREOST(IR0S) = IR0
    IREOTS(IR0) = IR0S
  end do
end do
NR0OB = IR0-NACOB
!write(u6,*) ' IR0 ',IR0

! RAS 4 orbitals

IR4 = NACOB+NR0OB
do ISM=1,NSMOB
  if (ISM == 1) then
    IBSM = 1
  else
    IBSM = IBSM+NTOOBS(ISM-1)
  end if
  NPREVS = NINOBS(ISM)+NR0OBS(ISM)+NACOBS(ISM)
  do ITP=1,NR4TP
    do IIR4=1,NR4OBS(ISM,ITP)
      ! Type ordered index
      IR4 = IR4+1
      ! Symmetry ordered index
      IR4S = IBSM+NPREVS-1+IIR4
      ISFTO(IR4) = ISM
      !ITFSO(IR4S) = 5
      IREOST(IR4S) = IR4
      IREOTS(IR4) = IR4S
    end do
  end do
end do
NR4OB = IR4-NACOB-NR0OB
!write(u6,*) ' IR4 ',IR4

! Inactive orbitals

IIN = NACOB+NR0OB+NR4OB
do ISM=1,NSMOB
  if (ISM == 1) then
    IBSM = 1
  else
    IBSM = IBSM+NTOOBS(ISM-1)
  end if
  NPREVS = 0
  do IIIN=1,NINOBS(ISM)
    ! Type ordered index
    IIN = IIN+1
    ! Symmetry ordered index
    IINS = IBSM+NPREVS-1+IIIN
    ISFTO(IIN) = ISM
    !ITFSO(IINS) = 6
    IREOST(IINS) = IIN
    IREOTS(IIN) = IINS
  end do
end do
NINOB = IIN-NACOB-NR0OB-NR4OB
!write(u6,*) ' IIN ',IIN

! Deleted orbitals

IDE = NACOB+NR0OB+NR4OB+NINOB
do ISM=1,NSMOB
  if (ISM == 1) then
    IBSM = 1
  else
    IBSM = IBSM+NTOOBS(ISM-1)
  end if
  IR4 = 0
  do ITP=1,NR4TP
    IR4 = IR4+NR4OBS(ISM,ITP)
  end do

  NPREVS = NINOBS(ISM)+NR0OBS(ISM)+NACOBS(ISM)+IR4
  do IIDE=1,NDEOBS(ISM)
    ! Type ordered index
    IDE = IDE+1
    ! Symmetry ordered index
    IDES = IBSM+NPREVS-1+IIDE
    ISFTO(IDE) = ISM
    !ITFSO(IDES) = 7
    IREOST(IDES) = IDE
    IREOTS(IDE) = IDES
  end do
end do
!write(u6,*) ' IDE ',IDE

IOFF = 1
do ISM=1,NSMOB
  IBSO(ISM) = IOFF
  IOFF = IOFF+NTOOBS(ISM)
end do

! ==================
! NTSOB,IBTSOB,ITSOB
! ==================

IOFF = 1
do I123=1,3
  do ISM=1,NSMOB
    NTSOB(I123,ISM) = NRSOBS(ISM,I123)
    IBTSOB(I123,ISM) = IOFF
    ITSOB(IOFF:IOFF+NRSOBS(ISM,I123)-1) = [(i,i=IOFF,IOFF+NRSOBS(ISM,I123)-1)]
    IOFF = IOFF+NRSOBS(ISM,I123)
  end do
end do

! ===================
! NOBPTS NOBPT IOBPTS
! ===================

! Loop over types in input order
call iCopy(NR4TP+6,[0],0,NOBPT,1)
LORB = 0   ! dummy initialize
IOTYPE = 0 ! dummy initialize
do ISMOB=1,NSMOB
  LSMOB = 0
  do ITYPE=1,NR4TP+6
    if (ITYPE == 1) then
      ! Inactive (frozen in normal notation)
      LORB = NINOBS(ISMOB)
      IOTYPE = 5+NR4TP
    else if (ITYPE == 2) then
      ! RAS0 (inactive in normal notation)
      LORB = NR0OBS(ISMOB)
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
    else if ((ITYPE >= 6) .and. (ITYPE <= 6+NR4TP-1)) then
      ! RAS4
      LORB = NR4OBS(ISMOB,ITYPE-5)
      IOTYPE = ITYPE-1
    else if (ITYPE == 6+NR4TP) then
      ! deleted orbitals
      LORB = NDEOBS(ISMOB)
      IOTYPE = ITYPE
    end if
    IOBPTS(IOTYPE,ISMOB) = LSMOB+1
    NOBPTS(IOTYPE,ISMOB) = LORB
    NOBPT(IOTYPE) = NOBPT(IOTYPE)+LORB
    LSMOB = LSMOB+LORB
  end do
end do

! ======
! ISMFSO
! ======

IORB = 0
do ISM=1,NSMOB
  do IOB=1,NTOOBS(ISM)
    IORB = IORB+1
    ISMFSO(IORB) = ISM
  end do
end do

#ifdef _DEBUGPRINT_
NTOOB = IDE
write(u6,*) ' =================='
write(u6,*) ' Output from ORBORD'
write(u6,*) ' =================='
write(u6,*) ' Symmetry of orbitals, type ordered'
call IWRTMA(ISFTO,1,NTOOB,1,NTOOB)
write(u6,*) ' Symmetry => type reordering array'
call IWRTMA(IREOST,1,NTOOB,1,NTOOB)
write(u6,*) ' Type => symmetry reordering array'
call IWRTMA(IREOTS,1,NTOOB,1,NTOOB)
write(u6,*) ' IBSO array'
call IWRTMA(IBSO,1,NSMOB,1,NSMOB)

write(u6,*) ' NTSOB array :'
call IWRTMA(NTSOB,3,NSMOB,3,NSMOB)
write(u6,*) ' IBTSOB array'
call IWRTMA(IBTSOB,3,NSMOB,3,NSMOB)
write(u6,*) ' ITSOB'
call IWRTMA(ITSOB,1,NACOB,1,NACOB)

write(u6,*) ' NOBPTS'
call IWRTMA(NOBPTS,6+NR4TP,NSMOB,6+MXPR4T,MXPOBS)
write(u6,*) ' NOBPT'
call IWRTMA(NOBPTS,6+NR4TP,1,6+MXPR4T,1)
write(u6,*) ' IOBPTS'
call IWRTMA(IOBPTS,6+NR4TP,NSMOB,6+MXPR4T,MXPOBS)

write(u6,*) ' ISFTO array :'
call IWRTMA(ISFTO,1,NTOOB,1,NTOOB)
!write(u6,*) ' ITFSO array :'
!call IWRTMA(ITFSO,1,NTOOB,1,NTOOB)

write(u6,*) ' ISMFSO array :'
call IWRTMA(ISMFSO,1,NTOOB,1,NTOOB)
write(u6,*) ' ITPFTO array :'
call IWRTMA(ITPFTO,1,NTOOB,1,NTOOB)
#endif

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer_array(NOCOBS)
  call Unused_integer_array(ITFSO)
end if

end subroutine ORBORD
