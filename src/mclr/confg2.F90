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
subroutine CONFG2(NORB1,NORB2,NORB3,NEL1MN,NEL3MX,MINOP,MAXOP,IREFSM,NEL,ICONF,NCNFTP)
! Generate array,ICONF,giving occupation of each configuration
! for CI space of reference symmetry IREFSM.
!
! Jeppe Olsen April 1989
!             August 1990 : Improved handling of large RAS 3 space
!
! Turbo configuration generator
! Iconf is ordered so all configurations of the same type are
! consecutively stored.
! ICONF is written so closed orbitals are given first and then single
! occupied orbitals

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NORB1, NORB2, NORB3, NEL1MN, NEL3MX, MINOP, MAXOP, IREFSM, NEL, NCNFTP(*)
integer(kind=iwp), intent(_OUT_) :: ICONF(*)
integer(kind=iwp) :: ICFREE, ICL, ICL1, IEL1, IEL1C, IEL3, IEL3C, IFRSTC, IFRSTO, IFSTR3, IIICHK, IOP, IORB, IORB1F, IORB1L, &
                     IORB2F, IORB2L, IORB3F, IORB3L, IPLACE, IPRORB, IR3CHK, ISYM, JCONF, K, KEL, KORB, MINCL1, MXMPTY, NCL, &
                     NEWORB, NOP, NORB
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: I, IBAS, IOC, IOPEN, ITYPE, LICONF
#endif
logical(kind=iwp) :: Loop700, Loop800, Skip700, Skip800, Test
integer(kind=iwp), allocatable :: IICL(:), IIOC(:), IIOP(:)
integer(kind=iwp), external :: ISYMST_MCLR

#ifndef _DEBUGPRINT_
#include "macros.fh"
unused_var(NCNFTP(1))
#endif

IORB1F = 1
IORB1L = IORB1F+NORB1-1

IORB2F = IORB1L+1
IORB2L = IORB2F+NORB2-1

IORB3F = IORB2L+1
IORB3L = IORB3F+NORB3-1

NORB = NORB1+NORB2+NORB3

call mma_allocate(IIOC,NORB,Label='IIOC')
call mma_allocate(IIOP,NORB,Label='IIOP')
call mma_allocate(IICL,NORB,Label='IICL')

! Loop over types of configurations

JCONF = 0
ICFREE = 1
! Min number of doubly occupied orbitals in RAS 1
MINCL1 = max(0,NEL1MN-NORB1)
#ifdef _DEBUGPRINT_
write(u6,*) ' Min number of doubly occupied orbitals in RAS 1',MINCL1
#endif
outer: do NOP=MINOP,MAXOP,2
  NCL = (NEL-NOP)/2
# ifdef _DEBUGPRINT_
  ITYPE = NOP-MINOP+1
  write(u6,*) ' NOP NCL ITYPE',NOP,NCL,ITYPE
# endif

  ! first combination of double occupied orbitals
  IIOC(1:NCL) = 2
  IIOC(NCL+1:NORB) = 0
  IICL(1:NCL) = [(ICL,ICL=1,NCL)]
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
        else if ((IPLACE /= NCL) .or. (NEWORB < NORB)) then

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
    do
      IORB = IORB-1
      if (IIOC(IORB) == 2) then
        MXMPTY = IORB-1
        if (IORB == 1) exit
      else
        exit
      end if
    end do
#   ifdef _DEBUGPRINT_
    write(u6,*) ' Next inactive configuration'
    call IWRTMA(IICL,1,NCL,1,NCL)
#   endif

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

          ! PAM 2013: Searching for next orbital with IIOC=0:
          do
            Test = NEWORB <= MXMPTY
            if (Test) Test = IIOC(NEWORB) /= 0
            if (.not. Test) exit
            NEWORB = NEWORB+1
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
              if ((IIOC(NEWORB) == 0) .or. (NEWORB >= MXMPTY)) exit
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
        ! RAS  CONSTRAINTS
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
        ! Faster routine for RAS 3, added august 1990
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
        ISYM = ISYMST_MCLR(IIOP,NOP)
        if (ISYM == IREFSM) then
#         ifdef _DEBUGPRINT_
          write(u6,1120) (IIOC(I),I=1,NORB)
#         endif
          JCONF = JCONF+1

          ICONF(ICFREE:ICFREE+NCL-1) = IICL(1:NCL)
          ICONF(ICFREE+NCL:ICFREE+NCL+NOP-1) = IIOP(1:NOP)
          ICFREE = ICFREE+NOP+NCL
        end if
      end if

      ! LOOP OVER active configurations end

      if ((NOP == 0) .and. (NCL == 0)) exit outer
      if (NOP == 0) exit
    end do
  end do inner
end do outer

call mma_deallocate(IIOC)
call mma_deallocate(IIOP)
call mma_deallocate(IICL)

#ifdef _DEBUGPRINT_
write(u6,'(/A,I3)') '  Configurations of symmetry ',IREFSM
write(u6,*) ' ================================='
IBAS = 0
do IOPEN=MINOP,MAXOP
  ITYPE = IOPEN-MINOP+1
  ICL = (NEL-IOPEN)/2
  IOC = IOPEN+ICL
  LICONF = NCNFTP(ITYPE)
  write(u6,'(/A,2I3)') '  Type with number of closed and open orbitals ',ICL,IOPEN
  write(u6,'(A,I7)') '  Number of configurations of this type',LICONF
  do JCONF=1,LICONF
    write(u6,'(3X,20I3)') (ICONF(IBAS+IORB),IORB=1,IOC)
    IBAS = IBAS+IOC
  end do
end do

1120 format('0  configuration included ',15I3)
#endif

end subroutine CONFG2
