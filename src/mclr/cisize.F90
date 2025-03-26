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
! Copyright (C) 1984,1989-1993, Jeppe Olsen                            *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine CISIZE(NORB1,NORB2,NORB3,NEL1MN,NEL3MX,NACTEL,MINOP,MAXOP,MXPCNT,MXPCSM,NCNATS,NCNASM,NDTASM,NCSASM,NDPCNT,NCPCNT,IICL, &
                  IIOP,IIOC)
! Number of configurations per per configuration type and symmetry
!
! Jeppe Olsen
!        August 1990 : Improved handling of large RAS 3 space
!        Winter 1991 : Modified for LUCIA

implicit none
integer NORB1, NORB2, NORB3, NEL1MN, NEL3MX, NACTEL, MINOP, MAXOP, MXPCNT, MXPCSM
! Output
integer NCNATS(MXPCNT,*), NCNASM(*), NDTASM(*), NCSASM(*)
! Input
integer NDPCNT(*), NCPCNT(*)
! Scratch
integer IICL(*), IIOP(*), IIOC(NORB1+NORB2+NORB3)
! Local variables
logical Test
integer ILOOP, ILOOP2, NCNF, NORBT, IORB1F, IORB1L, IORB2F, IORB2L, IORB3F, IORB3L, NORB, MINCL1, NOP, ITYPE, NCL, ICL, IFRSTC, &
        IORB, IPLACE, IPRORB, NEWORB, IEL1C, IEL3C, ICL1, IIICHK, MXMPTY, IOP, IFRSTO, IEL1, IEL3, IR3CHK, IFSTR3, K, KEL, KORB, &
        ISYM, NTYP, ISYMCN_MCLR
#ifdef _DEBUGPRINT_
integer I, ICSM
#endif

ILOOP = 0
ILOOP2 = 0
NCNF = 0
NORBT = NORB1+NORB2+NORB3

call iCOPY(MXPCSM*MXPCNT,[0],0,NCNATS,1)
call iCOPY(MXPCSM,[0],0,NCSASM,1)
call iCOPY(MXPCSM,[0],0,NDTASM,1)
call iCOPY(MXPCSM,[0],0,NCNASM,1)

IORB1F = 1
IORB1L = IORB1F+NORB1-1

IORB2F = IORB1L+1
IORB2L = IORB2F+NORB2-1

IORB3F = IORB2L+1
IORB3L = IORB3F+NORB3-1

NORB = NORB1+NORB2+NORB3
! Min number of doubly occupied orbitals in RAS 1
MINCL1 = max(0,NEL1MN-NORB1)
#ifdef _DEBUGPRINT_
write(6,*) ' Min number of doubly occupied orbitals in RAS 1',MINCL1
#endif
do NOP=MINOP,MAXOP,2
  ITYPE = NOP-MINOP+1
  NCL = (NACTEL-NOP)/2
# ifdef _DEBUGPRINT_
  write(6,*) ' NOP NCL ITYPE',NOP,NCL,ITYPE
# endif
  ! first combination of double occupied orbitals
  call iCOPY(NORB,[0],0,IIOC,1)
  do ICL=1,NCL
    IICL(ICL) = ICL
    IIOC(ICL) = 2
  end do
  IFRSTC = 1
  ! Loop over double occupied orbital configurations
2000 continue

  ! next double occupied configuration
  if ((IFRSTC == 1) .or. (NCL == 0)) goto 801
  !if ((IFRSTC == 0) .and. (NCL /= 0)) then

  do IORB=1,NORB
    if (IIOC(IORB) == 1) IIOC(IORB) = 0
  end do

  IPLACE = 0
800 continue
  IPLACE = IPLACE+1

  IPRORB = IICL(IPLACE)
  IIOC(IPRORB) = 0
  NEWORB = IPRORB+1
  if (((IPLACE < NCL) .and. (NEWORB < IICL(IPLACE+1))) .or. (IPLACE == NCL) .and. (NEWORB <= NORB)) then
    IICL(IPLACE) = NEWORB
    IIOC(NEWORB) = 2
  else if (.not. ((IPLACE == NCL) .and. (NEWORB >= NORB))) then
    if (IPLACE == 1) then
      IICL(1) = 1
      IIOC(1) = 2
    else
      IICL(IPLACE) = IICL(IPLACE-1)+1
      IIOC(IICL(IPLACE)) = 2
    end if
    goto 800
  else
    ! No more inactive configurations
    goto 2001
  end if
!OLD end if
801 continue
  IFRSTC = 0
# ifdef _DEBUGPRINT_
  write(6,*) ' Next inactive configuration'
  call IWRTMA(IICL,1,NCL,1,NCL)
# endif
  ! CHECK RAS1 and RAS 3
  IEL1C = 0
  IEL3C = 0
  ICL1 = 0
  do ICL=1,NCL
    IORB = IICL(ICL)
    if ((IORB1F <= IORB) .and. (IORB <= IORB1L)) then
      IEL1C = IEL1C+2
      ICL1 = ICL1+1
    else if ((IORB3F <= IORB) .and. (IORB <= IORB3L)) then
      IEL3C = IEL3C+2
    end if
  end do
  IIICHK = 1
  if ((ICL1 < MINCL1) .and. (IIICHK == 1)) then
    ! Next higher combination with a higher number of inactive orbitals
    do ICL=1,ICL1+1
      IIOC(IICL(ICL)) = 0
      IICL(ICL) = ICL
      IIOC(ICL) = 2
    end do
    IPLACE = ICL1+1
    if (IPLACE >= NCL) goto 2001
    goto 800
  end if
  if (IEL3C > NEL3MX) goto 2000
  ! Highest orbital not occupied
  MXMPTY = NORB
  IORB = NORB+1
  ! begin while
12 continue
  IORB = IORB-1
  if (IIOC(IORB) == 2) then
    MXMPTY = IORB-1
    if (IORB /= 1) goto 12
  end if
  ! End while

  ! first active configuration
  IORB = 0
  IOP = 0
  do IORB=1,NORB
    if (IIOC(IORB) == 0) then
      IOP = IOP+1
      if (IOP > NOP) goto 31
      IIOC(IORB) = 1
      IIOP(IOP) = IORB
    end if
  end do
31 continue
  IFRSTO = 1

  ! Next open shell configuration
1000 continue
  if ((IFRSTO == 1) .or. (NOP == 0)) goto 701
  !OLD if ((IFRSTO == 0) .and. (NOP /= 0)) then
  IPLACE = 0
700 continue
  IPLACE = IPLACE+1
  IPRORB = IIOP(IPLACE)
  NEWORB = IPRORB+1
  IIOC(IPRORB) = 0
690 continue
  if ((NEWORB <= MXMPTY) .and. (IIOC(min(NORBT,NEWORB)) /= 0)) then
    NEWORB = NEWORB+1
    if (NEWORB <= MXMPTY) goto 690
  end if
  !691 End of loop
  Test = IPLACE < NOP
  if (Test) Test = NEWORB < IIOP(IPLACE+1)
  if (Test .or. (IPLACE == NOP) .and. (NEWORB <= MXMPTY)) then
    IIOP(IPLACE) = NEWORB
    IIOC(NEWORB) = 1
  else if (IPLACE /= NOP) then
    if (IPLACE == 1) then
      NEWORB = 1-1
    else
      NEWORB = IIOP(IPLACE-1)
    end if
671 continue
    NEWORB = NEWORB+1
    if ((IIOC(NEWORB) /= 0) .and. (NEWORB <= MXMPTY)) goto 671
    IIOP(IPLACE) = NEWORB
    IIOC(NEWORB) = 1
    goto 700
  else
    ! No more active configurations, so
    if (NCL /= 0) goto 2000
    if (NCL == 0) goto 5001
  end if
!end if
701 continue
  IFRSTO = 0

# ifdef _DEBUGPRINT_
  write(6,*) ' Next active configuration'
  call IWRTMA(IIOP,1,NOP,1,NOP)
# endif
  ! RAS CONSTRAINTS
  IEL1 = IEL1C
  IEL3 = IEL3C
  ! CHECK RAS1 and RAS3
  do IOP=1,NOP
    IORB = IIOP(IOP)
    if ((IORB1F <= IORB) .and. (IORB <= IORB1L)) then
      IEL1 = IEL1+1
    else if ((IORB3F <= IORB) .and. (IORB <= IORB3L)) then
      IEL3 = IEL3+1
    end if
  end do
  IR3CHK = 1
  if ((IEL3 > NEL3MX) .and. (IR3CHK == 1)) then
    ! Number of electrons in substring
    IFSTR3 = 0
    do IOP=1,NOP
      if (IIOP(IOP) >= IORB3F) then
        IFSTR3 = IOP
        goto 5608
      end if
    end do
5608 continue
    if (IFSTR3 /= NOP) then

      ! Lowest possible string with NOP electrons
      do K=1,IFSTR3
        IIOC(IIOP(K)) = 0
      end do

      KEL = 0
      KORB = 0
5630  continue
      KORB = KORB+1
      if (IIOC(KORB) /= 2) then
        KEL = KEL+1
        IIOC(KORB) = 1
        IIOP(KEL) = KORB
      end if
      if (KEL /= IFSTR3) goto 5630
      IPLACE = IFSTR3
      goto 700
    end if
  end if
  if ((IEL1 < NEL1MN) .or. (IEL3 > NEL3MX)) goto 999
  ! Spatial symmetry
  ISYM = ISYMCN_MCLR(IICL,IIOP,NCL,NOP)

# ifdef _DEBUGPRINT_
  write(6,*) ' ISYM : ',ISYM
  write(6,1120) (IIOC(I),I=1,NORB)
# endif
  NCNF = NCNF+1

  NCNASM(ISYM) = NCNASM(ISYM)+1
# ifdef _DEBUGPRINT_
  write(6,1311) NCNF,(IIOC(I),I=1,NORB)
# endif

  NCNATS(ITYPE,ISYM) = NCNATS(ITYPE,ISYM)+1
# ifdef _DEBUGPRINT_
  write(6,3111) NCNF,ITYPE
# endif

! LOOP OVER CONFIGURATIONS, end

999 continue
  ILOOP = ILOOP+1
  ILOOP2 = ILOOP2+1
  if (ILOOP2 == 10000000) then
    write(6,*) ' 10 million configurations generated'
    ILOOP2 = 0
  end if

  if ((NOP == 0) .and. (NCL == 0)) goto 5001
  if (NOP == 0) goto 2000
  goto 1000
2001 continue
end do
5001 continue

#ifdef _DEBUGPRINT_
write(6,'(A,I8)') '  Total number of configurations generated ',NCNF
#endif
! ==============================
! Total number of CSF's and SD's
! ==============================
NTYP = MAXOP-MINOP+1
do ISYM=1,MXPCSM
  do ITYPE=1,NTYP
    NDTASM(ISYM) = NDTASM(ISYM)+NDPCNT(ITYPE)*NCNATS(ITYPE,ISYM)
    NCSASM(ISYM) = NCSASM(ISYM)+NCPCNT(ITYPE)*NCNATS(ITYPE,ISYM)
  end do
end do

#ifdef _DEBUGPRINT_
write(6,*)
write(6,'(/A)') ' Information about actual configurations'
write(6,'(A)') ' ========================================'
write(6,'(/A)') '    Symmetry     Configurations     CSFs     Combinations'
write(6,'(A)') '  =============  ============== ============ ============'
do ICSM=1,MXPCSM
  if (NCNASM(ICSM) /= 0) write(6,'(4X,I3,4X,6X,I8,6X,I8,6X,I9)') ICSM,NCNASM(ICSM),NCSASM(ICSM),NDTASM(ICSM)
end do

return

1120 format('0  configuration included ',15I3,('                         ',15I3))
1311 format('  configuration ',I3,20I2,/,(1X,18X,20I2))
3111 format('0  CONFIGURATION..',I3,' IS TYPE..',I3)
#endif

end subroutine CISIZE
