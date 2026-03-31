!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine SYG2SGU(IMODE,SGS,CIS,LSYM,ICNFTAB,ISPNTAB,CIOLD,CINEW)
! SGS       : Data that define a Split Graph
! qCIS : Data that define a CI array structure
! IMODE=0 transforms a Symmetric Group CI array to SGUGA
! IMODE=1 transforms a Split GUGA CI array to Symm Group
! ...Configuration and Spin Coupling tables, fill this in later.
! CIOLD and CINEW are obvious.

use rassi_aux, only: ipglob
use gugx, only: CIStruct, SGStruct
use Molcas, only: MxLev
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IMODE, ICNFTAB(*), ISPNTAB(*)
type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(in) :: CIS
integer(kind=iwp), intent(inout) :: LSYM
real(kind=wp), intent(in) :: CIOLD(*)
real(kind=wp), intent(out) :: CINEW(*)
integer(kind=iwp), parameter :: IFUP2CS(0:1) = [2,1], MXCPI = 15, NBUFFER = 600
integer(kind=iwp) :: I, ICASE(400), ICNF, ICNUM(NBUFFER), ICPL, ICSPLT, ICSYMG, IEL, IEL1, IEL2, IFORM, IFUP, IOCC, IORB, IREST, &
                     IWLKPOS, IWORD, IWRD, KCNF, KCNFINF, KCPL, KGSLIM, KGSORB, KSPNINF, KWALK(NBUFFER), LEV, MAXOP, MINOP, &
                     MIPWLK, MXWLK, NACTEL, NAPART, NCLSD, NCNF, NCONF, NCPL, NCSYMG, NHEAD, NLEV, NOCC, NODD, NOPEN, NSYM, NWALK, &
                     NWLKLST, NWRD
real(kind=wp) :: PHASE(NBUFFER), PHS
integer(kind=iwp), allocatable :: MWS2W(:), OrbArr(:)

! Dereference CIS and SGS       for some data:
NCONF = CIS%NCSF(LSYM)
NWALK = CIS%nWalk
call mma_allocate(MWS2W,NWALK,Label='MWS2W')
NSYM = ICNFTAB(7)
call MSTOW(SGS,CIS,MWS2W,NSYM)
! MWS2W is a table which gives the upper or lower walk
! index as function of the MAW sum.

! Inspect the top row of the DRT to find NACTEL and spin:
NLEV = SGS%DRT(1,1)
if (NLEV > MXLEV) then
  write(u6,*) ' SYG2SGU: error: number of levels exceeds MXLEV'
  write(u6,'(1X,2(A,I4))') ' NLEV = ',NLEV,' MXLEV = ',MXLEV
  call AbEnd()
end if
NACTEL = SGS%DRT(1,2)

! Now a good bound on MINOP, the minimum number of open
! shells, would be MLTPLC-1. This is the best bound, and it
! does not depend on any assumed Ms.

! A buffer of packed walks is used:
NWLKLST = 0
IWLKPOS = 1
! Nr of integers used to store each total walk:
MIPWLK = 1+(NLEV-1)/MXCPI
! Max nr of Split-Graph CSF''s, stored as walks in the
! buffer KWALK
MXWLK = NBUFFER/MIPWLK
! Pick up various data from conf and spin tables
! Unbutton Configuration table:

MINOP = ICNFTAB(5)
MAXOP = ICNFTAB(6)
NSYM = ICNFTAB(7)
if (LSYM /= ICNFTAB(8)) call abend()
LSYM = ICNFTAB(8)
NAPART = ICNFTAB(9)
IFORM = ICNFTAB(10)
!PAM07 Statement to fool intel 10.1 compiler do do the right thing:
if (MINOP > MAXOP) write(u6,*) MINOP,MAXOP

NHEAD = 10
KGSORB = NHEAD+1
KGSLIM = KGSORB+(NSYM+1)*(NAPART+1)
KCNFINF = KGSLIM+2*NAPART
! Unbutton Spin table:
KSPNINF = 9
call mma_allocate(ORBARR,NACTEL,Label='OrbArr')

! Loop over nr of open shells
! NCSYMG=Nr of Symmetric-Group CSF''s treated so far.
NCSYMG = 0
do NOPEN=MINOP,MAXOP
  NCLSD = (NACTEL-NOPEN)/2
  if (NCLSD < 0) cycle
  if (2*NCLSD+NOPEN /= NACTEL) cycle
  NOCC = NCLSD+NOPEN
  if (NOCC > NLEV) cycle
  NCNF = ICNFTAB(KCNFINF+3*(LSYM-1+NSYM*(NOPEN-MINOP)))
  KCNF = ICNFTAB(KCNFINF+3*(LSYM-1+NSYM*(NOPEN-MINOP))+1)
  NWRD = ICNFTAB(KCNFINF+3*(LSYM-1+NSYM*(NOPEN-MINOP))+2)
  NCPL = ISPNTAB(KSPNINF+6*(NOPEN-MINOP)+1)
  KCPL = ISPNTAB(KSPNINF+6*(NOPEN-MINOP)+3)

  ! Here follows four if-clauses on the cases of IFORM=1..4.
  ! The contents in all four cases is almost the same, but it is
  ! *DELIBERATE* that it has not been rewritten as one single piece of code
  ! with a few IF statement inside. It seems that some compilers do not
  ! optimize the code correctly then.
  ! Also, a number of IF-clauses have been replaced by arithmetic computation
  ! for the same reason. Pardon the clumsy result...
  IWORD = 0 ! dummy initialize
  if (IFORM == 1) then
    !*******************************************************************
    ! The IFORM=1 case:
    ! Long loop over configurations
    do ICNF=1,NCNF
      do IEL=1,NOCC
        ORBARR(IEL) = ICNFTAB(KCNF-1+IEL+NWRD*(ICNF-1))
      end do
      ! Loop over spin couplings
      do ICPL=1,NCPL
        do I=1,NLEV
          ICASE(I) = 0
        end do
        do I=1,NCLSD
          IORB = ORBARR(I)
          ICASE(IORB) = 3
        end do
        do I=1,NOPEN
          IORB = ORBARR(NCLSD+I)
          IFUP = ISPNTAB(KCPL-1+I+NOPEN*(ICPL-1))
          !ICASE(IORB) = 1
          !if (IFUP /= 1) ICASE(IORB) = 2
          ICASE(IORB) = IFUP2CS(IFUP)

        end do
        ! A phase factor will be induced by the reordering.
        PHS = One
        NODD = 0
        do LEV=1,NLEV
          I = ICASE(LEV)
          !if (I == 1) NODD = 1-NODD
          !if (I == 2) NODD = 1-NODD
          NODD = ((2*NODD-1)*I*(I-3))/2+NODD
          !if (NODD == 1) then
          !  IF(I == 2) PHS = -PHS
          !  IF(I == 3) PHS = -PHS
          !end if
          PHS = real((3+NODD*I*(I-1)*(2*I-7))/3,kind=wp)*PHS
        end do
        ! Pack the walk and add it to the list.
        call PKWLK(NLEV,MIPWLK,1,KWALK(IWLKPOS),ICASE)
        IWLKPOS = IWLKPOS+MIPWLK
        NWLKLST = NWLKLST+1
        PHASE(NWLKLST) = PHS
        if (NWLKLST == MXWLK) then
          ! Translate the packed walks to a list of CSF ID-numbers
          call W2SGORD(SGS,CIS,MWS2W,NWLKLST,KWALK,ICNUM)
          IWLKPOS = 1
          ! Loop over this list.
          ICSYMG = NCSYMG
          if (IMODE == 0) then
            do I=1,NWLKLST
              PHS = PHASE(I)
              ICSPLT = ICNUM(I)
              ICSYMG = NCSYMG+I
              CINEW(ICSPLT) = CIOLD(ICSYMG)*PHS
            end do
          else
            do I=1,NWLKLST
              PHS = PHASE(I)
              ICSPLT = ICNUM(I)
              ICSYMG = NCSYMG+I
              CINEW(ICSYMG) = CIOLD(ICSPLT)*PHS
            end do
          end if
          NCSYMG = ICSYMG
          NWLKLST = 0
        end if

        ! End of loop over spin couplings
      end do
      ! End of loop over configurations
    end do
  else if (IFORM == 2) then
    !*******************************************************************
    ! The IFORM=2 case:
    ! Long loop over configurations
    do ICNF=1,NCNF
      IEL2 = 0
      IEL1 = NCLSD
      do IORB=1,NLEV
        IOCC = ICNFTAB(KCNF-1+IORB+NWRD*(ICNF-1))
        if (IOCC == 1) then
          IEL1 = IEL1+1
          ORBARR(IEL1) = IORB
        else
          IEL2 = IEL2+1
          ORBARR(IEL2) = IORB
        end if
      end do
      ! Loop over spin couplings
      do ICPL=1,NCPL
        do I=1,NLEV
          ICASE(I) = 0
        end do
        do I=1,NCLSD
          IORB = ORBARR(I)
          ICASE(IORB) = 3
        end do
        do I=1,NOPEN
          IORB = ORBARR(NCLSD+I)
          IFUP = ISPNTAB(KCPL-1+I+NOPEN*(ICPL-1))
          !ICASE(IORB) = 1
          !if (IFUP /= 1) ICASE(IORB) = 2
          ICASE(IORB) = IFUP2CS(IFUP)
        end do
        ! A phase factor will be induced by the reordering.
        PHS = One
        NODD = 0
        do LEV=1,NLEV
          I = ICASE(LEV)
          !if (I == 1) NODD = 1-NODD
          !if (I == 2) NODD = 1-NODD
          NODD = ((2*NODD-1)*I*(I-3))/2+NODD
          !if (NODD == 1) then
          !  if (I == 2) PHS = -PHS
          !  if (I == 3) PHS = -PHS
          !end if
          PHS = real((3+NODD*I*(I-1)*(2*I-7))/3,kind=wp)*PHS
        end do
        ! Pack the walk and add it to the list.
        call PKWLK(NLEV,MIPWLK,1,KWALK(IWLKPOS),ICASE)
        IWLKPOS = IWLKPOS+MIPWLK
        NWLKLST = NWLKLST+1
        PHASE(NWLKLST) = PHS
        if (NWLKLST == MXWLK) then
          ! Translate the packed walks to a list of CSF ID-numbers
          call W2SGORD(SGS,CIS,MWS2W,NWLKLST,KWALK,ICNUM)
          IWLKPOS = 1
          ! Loop over this list.
          ICSYMG = NCSYMG
          if (IMODE == 0) then
            do I=1,NWLKLST
              PHS = PHASE(I)
              ICSPLT = ICNUM(I)
              ICSYMG = NCSYMG+I
              CINEW(ICSPLT) = CIOLD(ICSYMG)*PHS
            end do
          else
            do I=1,NWLKLST
              PHS = PHASE(I)
              ICSPLT = ICNUM(I)
              ICSYMG = NCSYMG+I
              CINEW(ICSYMG) = CIOLD(ICSPLT)*PHS
            end do
          end if
          NCSYMG = ICSYMG
          NWLKLST = 0
        end if

        ! End of loop over spin couplings
      end do
      ! End of loop over configurations
    end do
  else if (IFORM == 3) then
    !*******************************************************************
    ! The IFORM=3 case:
    ! Long loop over configurations
    do ICNF=1,NCNF
      do IEL=1,NOCC
        IWRD = (3+IEL)/4
        IREST = (3+IEL)-4*IWRD
        if (IREST == 0) IWORD = ICNFTAB(KCNF-1+IWRD+NWRD*(ICNF-1))
        IORB = mod(IWORD,256)
        IWORD = IWORD/256
        ORBARR(IEL) = IORB
      end do
      ! Loop over spin couplings
      do ICPL=1,NCPL
        do I=1,NLEV
          ICASE(I) = 0
        end do
        do I=1,NCLSD
          IORB = ORBARR(I)
          ICASE(IORB) = 3
        end do
        do I=1,NOPEN
          IORB = ORBARR(NCLSD+I)
          IFUP = ISPNTAB(KCPL-1+I+NOPEN*(ICPL-1))
          !ICASE(IORB) = 1
          !if (IFUP /= 1) ICASE(IORB) = 2
          ICASE(IORB) = IFUP2CS(IFUP)
        end do
        ! A phase factor will be induced by the reordering.
        PHS = One
        NODD = 0
        do LEV=1,NLEV
          I = ICASE(LEV)
          !if (I == 1) NODD = 1-NODD
          !if (I == 2) NODD = 1-NODD
          NODD = ((2*NODD-1)*I*(I-3))/2+NODD
          !if (NODD == 1) then
          !  if (I == 2) PHS = -PHS
          !  if (I == 3) PHS = -PHS
          !end if
          PHS = real((3+NODD*I*(I-1)*(2*I-7))/3,kind=wp)*PHS
        end do
        ! Pack the walk and add it to the list.
        call PKWLK(NLEV,MIPWLK,1,KWALK(IWLKPOS),ICASE)
        IWLKPOS = IWLKPOS+MIPWLK
        NWLKLST = NWLKLST+1
        PHASE(NWLKLST) = PHS
        if (NWLKLST == MXWLK) then
          ! Translate the packed walks to a list of CSF ID-numbers
          call W2SGORD(SGS,CIS,MWS2W,NWLKLST,KWALK,ICNUM)
          IWLKPOS = 1
          ! Loop over this list.
          ICSYMG = NCSYMG
          if (IMODE == 0) then
            do I=1,NWLKLST
              PHS = PHASE(I)
              ICSPLT = ICNUM(I)
              ICSYMG = NCSYMG+I
              CINEW(ICSPLT) = CIOLD(ICSYMG)*PHS
            end do
          else
            do I=1,NWLKLST
              PHS = PHASE(I)
              ICSPLT = ICNUM(I)
              ICSYMG = NCSYMG+I
              CINEW(ICSYMG) = CIOLD(ICSPLT)*PHS
            end do
          end if
          NCSYMG = ICSYMG
          NWLKLST = 0
        end if

        ! End of loop over spin couplings
      end do
      ! End of loop over configurations
    end do
  else if (IFORM == 4) then
    !*******************************************************************
    ! The IFORM=4 case:
    ! Long loop over configurations
    do ICNF=1,NCNF
      IEL2 = 0
      IEL1 = NCLSD
      do IORB=1,NLEV
        IWRD = (IORB+14)/15
        IREST = IORB+14-15*IWRD
        if (IREST == 0) IWORD = ICNFTAB(KCNF-1+IWRD+NWRD*(ICNF-1))
        IOCC = mod(IWORD,4)
        IWORD = IWORD/4
        if (IOCC == 1) then
          IEL1 = IEL1+1
          ORBARR(IEL1) = IORB
        else
          IEL2 = IEL2+1
          ORBARR(IEL2) = IORB
        end if
      end do
      ! Loop over spin couplings
      do ICPL=1,NCPL
        do I=1,NLEV
          ICASE(I) = 0
        end do
        do I=1,NCLSD
          IORB = ORBARR(I)
          ICASE(IORB) = 3
        end do
        do I=1,NOPEN
          IORB = ORBARR(NCLSD+I)
          IFUP = ISPNTAB(KCPL-1+I+NOPEN*(ICPL-1))
          !ICASE(IORB) = 1
          !if (IFUP /= 1) ICASE(IORB) = 2
          ICASE(IORB) = IFUP2CS(IFUP)
        end do
        ! A phase factor will be induced by the reordering.
        PHS = One
        NODD = 0
        do LEV=1,NLEV
          I = ICASE(LEV)
          !if (I == 1) NODD = 1-NODD
          !if (I == 2) NODD = 1-NODD
          NODD = ((2*NODD-1)*I*(I-3))/2+NODD
          !if (NODD == 1) then
          !  if (I == 2) PHS = -PHS
          !  if (I == 3) PHS = -PHS
          !end if
          PHS = real((3+NODD*I*(I-1)*(2*I-7))/3,kind=wp)*PHS
        end do
        ! Pack the walk and add it to the list.
        call PKWLK(NLEV,MIPWLK,1,KWALK(IWLKPOS),ICASE)
        IWLKPOS = IWLKPOS+MIPWLK
        NWLKLST = NWLKLST+1
        PHASE(NWLKLST) = PHS
        if (NWLKLST == MXWLK) then
          ! Translate the packed walks to a list of CSF ID-numbers
          call W2SGORD(SGS,CIS,MWS2W,NWLKLST,KWALK,ICNUM)
          IWLKPOS = 1
          ! Loop over this list.
          ICSYMG = NCSYMG
          if (IMODE == 0) then
            do I=1,NWLKLST
              PHS = PHASE(I)
              ICSPLT = ICNUM(I)
              ICSYMG = NCSYMG+I
              CINEW(ICSPLT) = CIOLD(ICSYMG)*PHS
            end do
          else
            do I=1,NWLKLST
              PHS = PHASE(I)
              ICSPLT = ICNUM(I)
              ICSYMG = NCSYMG+I
              CINEW(ICSYMG) = CIOLD(ICSPLT)*PHS
            end do
          end if
          NCSYMG = ICSYMG
          NWLKLST = 0
        end if

        ! End of loop over spin couplings
      end do
      ! End of loop over configurations
    end do
  end if
  !*********************************************************************

  ! End of loop over NOPEN
end do

! As above, processing what remains in the KWALK buffer.
call W2SGORD(SGS,CIS,MWS2W,NWLKLST,KWALK,ICNUM)
IWLKPOS = 1
if (IMODE == 0) then
  do I=1,NWLKLST
    PHS = PHASE(I)
    ICSPLT = ICNUM(I)
    ICSYMG = NCSYMG+I
    CINEW(ICSPLT) = CIOLD(ICSYMG)*PHS
  end do
else
  do I=1,NWLKLST
    PHS = PHASE(I)
    ICSPLT = ICNUM(I)
    ICSYMG = NCSYMG+I
    CINEW(ICSYMG) = CIOLD(ICSPLT)*PHS
  end do
end if
NCSYMG = NCSYMG+NWLKLST
NWLKLST = 0

if (IPGLOB >= 5) then
  write(u6,*)
  write(u6,*) ' CI vector reordered in SYG2SGU'
  if (IMODE == 0) then
    write(u6,*) ' OLD CI vector in SYMMETRIC GROUP order'
    write(u6,'(10F12.8)') (CIOLD(I),I=1,NCONF)
    write(u6,*) ' NEW CI vector in SPLIT GRAPH UGA order'
    write(u6,'(10F12.8)') (CINEW(I),I=1,NCONF)
  else
    write(u6,*) ' OLD CI vector in SPLIT GRAPH UGA order'
    write(u6,'(10F12.8)') (CIOLD(I),I=1,NCONF)
    write(u6,*) ' NEW CI vector in SYMMETRIC GROUP order'
    write(u6,'(10F12.8)') (CINEW(I),I=1,NCONF)
  end if
end if

call mma_deallocate(OrbArr)
call mma_deallocate(MWS2W)

end subroutine SYG2SGU
