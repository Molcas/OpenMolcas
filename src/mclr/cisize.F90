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

use Definitions, only: u6

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
logical Loop700, Loop800, Skip700, Skip800

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
write(u6,*) ' Min number of doubly occupied orbitals in RAS 1',MINCL1
#endif
outer: do NOP=MINOP,MAXOP,2
  ITYPE = NOP-MINOP+1
  NCL = (NACTEL-NOP)/2
# ifdef _DEBUGPRINT_
  write(u6,*) ' NOP NCL ITYPE',NOP,NCL,ITYPE
# endif
  ! first combination of double occupied orbitals
  call iCOPY(NORB,[0],0,IIOC,1)
  do ICL=1,NCL
    IICL(ICL) = ICL
    IIOC(ICL) = 2
  end do
  IFRSTC = 1

  ! Loop over double occupied orbital configurations
  inner: do

    ! next double occupied configuration
    if ((IFRSTC == 1) .or. (NCL == 0)) then
      Skip800 = .true.
    else
      do IORB=1,NORB
        if (IIOC(IORB) == 1) IIOC(IORB) = 0
      end do
      IPLACE = 0
      Skip800 = .false.
    end if
    Loop800 = .true.
    do while (Loop800)
      Loop800 = .false.
      if (Skip800) then
        Skip800 = .false.
      else
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
          Loop800 = .true.
          cycle
        else
          ! No more inactive configurations
          exit inner
        end if
      end if
      IFRSTC = 0
#     ifdef _DEBUGPRINT_
      write(u6,*) ' Next inactive configuration'
      call IWRTMA(IICL,1,NCL,1,NCL)
#     endif
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
        if (IPLACE >= NCL) exit inner
        Loop800 = .true.
      end if
    end do
    if (IEL3C > NEL3MX) cycle inner
    ! Highest orbital not occupied
    MXMPTY = NORB
    IORB = NORB+1
    ! begin while
    do
      IORB = IORB-1
      if (IIOC(IORB) == 2) then
        MXMPTY = IORB-1
        if (IORB == 1) exit
      else
        exit
      end if
    end do

    ! first active configuration
    IORB = 0
    IOP = 0
    do IORB=1,NORB
      if (IIOC(IORB) == 0) then
        IOP = IOP+1
        if (IOP > NOP) exit
        IIOC(IORB) = 1
        IIOP(IOP) = IORB
      end if
    end do
    IFRSTO = 1

    ! Next open shell configuration
    do
      if ((IFRSTO == 1) .or. (NOP == 0)) then
        Skip700 = .true.
      else
        IPLACE = 0
        Skip700 = .false.
      end if
      Loop700 = .true.
      do while (Loop700)
        Loop700 = .false.
        if (Skip700) then
          Skip700 = .false.
        else
          IPLACE = IPLACE+1
          IPRORB = IIOP(IPLACE)
          NEWORB = IPRORB+1
          IIOC(IPRORB) = 0
          do
            if ((NEWORB <= MXMPTY) .and. (IIOC(min(NORBT,NEWORB)) /= 0)) then
              NEWORB = NEWORB+1
              if (NEWORB > MXMPTY) exit
            else
              exit
            end if
          end do
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
            do
              NEWORB = NEWORB+1
              if ((IIOC(NEWORB) == 0) .or. (NEWORB > MXMPTY)) exit
            end do
            IIOP(IPLACE) = NEWORB
            IIOC(NEWORB) = 1
            Loop700 = .true.
            cycle
          else
            ! No more active configurations, so
            if (NCL /= 0) cycle inner
            if (NCL == 0) exit outer
          end if
        end if
        IFRSTO = 0

#       ifdef _DEBUGPRINT_
        write(u6,*) ' Next active configuration'
        call IWRTMA(IIOP,1,NOP,1,NOP)
#       endif
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
              exit
            end if
          end do
          if (IFSTR3 /= NOP) then

            ! Lowest possible string with NOP electrons
            do K=1,IFSTR3
              IIOC(IIOP(K)) = 0
            end do

            KEL = 0
            KORB = 0
            do
              KORB = KORB+1
              if (IIOC(KORB) /= 2) then
                KEL = KEL+1
                IIOC(KORB) = 1
                IIOP(KEL) = KORB
              end if
              if (KEL == IFSTR3) exit
            end do
            IPLACE = IFSTR3
            Loop700 = .true.
          end if
        end if
      end do
      if ((IEL1 >= NEL1MN) .and. (IEL3 <= NEL3MX)) then
        ! Spatial symmetry
        ISYM = ISYMCN_MCLR(IICL,IIOP,NCL,NOP)

#       ifdef _DEBUGPRINT_
        write(u6,*) ' ISYM : ',ISYM
        write(u6,1120) (IIOC(I),I=1,NORB)
#       endif
        NCNF = NCNF+1

        NCNASM(ISYM) = NCNASM(ISYM)+1
#       ifdef _DEBUGPRINT_
        write(u6,1311) NCNF,(IIOC(I),I=1,NORB)
#       endif

        NCNATS(ITYPE,ISYM) = NCNATS(ITYPE,ISYM)+1
#       ifdef _DEBUGPRINT_
        write(u6,3111) NCNF,ITYPE
#       endif
      end if

      ! LOOP OVER CONFIGURATIONS, end

      ILOOP = ILOOP+1
      ILOOP2 = ILOOP2+1
      if (ILOOP2 == 10000000) then
        write(u6,*) ' 10 million configurations generated'
        ILOOP2 = 0
      end if

      if ((NOP == 0) .and. (NCL == 0)) exit outer
      if (NOP == 0) exit
    end do
  end do inner
end do outer

#ifdef _DEBUGPRINT_
write(u6,'(A,I8)') '  Total number of configurations generated ',NCNF
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
write(u6,*)
write(u6,'(/A)') ' Information about actual configurations'
write(u6,'(A)') ' ========================================'
write(u6,'(/A)') '    Symmetry     Configurations     CSFs     Combinations'
write(u6,'(A)') '  =============  ============== ============ ============'
do ICSM=1,MXPCSM
  if (NCNASM(ICSM) /= 0) write(u6,'(4X,I3,4X,6X,I8,6X,I8,6X,I9)') ICSM,NCNASM(ICSM),NCSASM(ICSM),NDTASM(ICSM)
end do

return

1120 format('0  configuration included ',15I3,('                         ',15I3))
1311 format('  configuration ',I3,20I2,/,(1X,18X,20I2))
3111 format('0  CONFIGURATION..',I3,' IS TYPE..',I3)
#endif

end subroutine CISIZE
