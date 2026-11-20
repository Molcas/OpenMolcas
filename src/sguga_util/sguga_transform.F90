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
! Copyright (C) 1989,1998, Per Ake Malmqvist                           *
!***********************************************************************

#include "intent.fh"

!> SGUGA: transformations of CI vectors of a split graph: the CI transformation
!> for orbital rotations (SG_CITRA) and the conversion between the symmetric
!> group and the split-graph CSF orders (SYG2SGU), with their walk tools.
module SGUGA_Transform

implicit none
private

public :: SG_CITRA, SYG2SGU

contains

!  PROGRAM RASSI        PER-AAKE MALMQVIST
!  SUBROUTINE SG_CITRA     IBM-3090 RELEASE 89 01 31
!  USE THE COEFFICIENTS FOR A SEQUENCE OF SINGLE-ORBITAL TRANSFOR-
!  MATION, TRA, TO TRANSFORM THE CI EXPANSION COEFFICIENTS
!  IN-PLACE TO A NEW NON-ON ORBITAL BASIS.
!  NEW VERSION 981122, using user define types SGS,CIS,XS.
!***********************************************************************
!  SG_CITRA
!
!> @brief
!>   Recompute a CI coefficient array to use another orbital basis
!> @author P. &Aring;. Malmqvist
!>
!> @details
!> For a given linear transformation of the orbitals, and a
!> CI array where the CSF basis is built from the old
!> orbitals, compute the CI array using the new orbitals
!> instead. The orbitals are transformed sequentially, and
!> for each active orbital, a call to ::SG_SSOTRA performs the
!> single-orbital transformation.
!>
!> @param[in]     WFTP Wave function Type Name
!> @param[in]     SGS Split Graph Structure user defined type
!> @param[in]     CIS CI Structure user define type
!> @param[in,out] EXS  Excitation operator Structure user defined type
!> @param[in]     LSM  Wave function Symmetry Label
!> @param[in]     TRA  Transformation Matrix
!> @param[in]     NCO  Number of Configuration Functions
!> @param[in,out] CI   CI Array
!***********************************************************************

!#ifdef _DEBUGPRINT_
subroutine SG_CITRA(WFTP,iState,LSM,NTRA,TRA,NCO,CI,NOSH,NISH,NASH)

  use Molcas, only: MxSym
  use Symmetry_Info, only: nIrrep
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero, One, Three, Half
  use Definitions, only: wp, iwp, u6

  implicit none
  character(len=8), intent(in) :: WFTP
  integer(kind=iwp), intent(in) :: iState, LSM, NCO, NTRA
  integer(kind=iwp), intent(in) :: NOSH(MxSym) ,NISH(MxSym) ,NASH(MxSym)
  real(kind=wp), intent(in) :: TRA(NTRA)
  real(kind=wp), intent(inout) :: CI(NCO)
  integer(kind=iwp) :: I, II, ISTA, ISYM, NA, NI, NO
  real(kind=wp) :: CKK, FAC
  real(kind=wp), allocatable :: TMP(:)
#ifdef _DEBUGPRINT_
  real(kind=wp), external :: ddot_
#endif

#ifdef _DEBUGPRINT_
  write(u6,*) ' Entering CITRA. norm=',ddot_(NCO,CI,1,CI,1)
#endif
  ! TRA contains square matrices, one per symmetry
  !  FIRST TRANSFORM THE INACTIVE ORBITALS:

  FAC = One
  ISTA = 1
  do ISYM=1,nIrrep
    NO = NOSH(ISYM)
    do I=1,NISH(ISYM)
      II = ISTA+(NO+1)*(I-1)
      CKK = TRA(II)
      FAC = FAC*CKK
    end do
    ISTA = ISTA+NO**2
  end do

  FAC = FAC**2
  CI(:) = FAC*CI(:)
  ! THEN THE ACTIVE ONES:
  Select Case(WFTP)

    Case ('HISPIN','CLOSED')
      ! The HISPIN case may be buggy and is not presently used.
      ISTA = 1
      do ISYM=1,nIrrep
        NI = NISH(ISYM)
        NA = NASH(ISYM)
        NO = NOSH(ISYM)
        do I=NI+1,NI+NA
          II = ISTA+(NO+1)*(I-1)
          CKK = TRA(II)
          FAC = FAC*CKK
        end do
        ISTA = ISTA+NO**2
      end do
      if (WFTP == 'CLOSED') FAC = FAC**2
      CI(:) = FAC*CI(:)

    Case ('EMPTY')

    Case ('GENERAL')
      ! The general case:
      call mma_allocate(TMP,NCO,Label='TMP')
      ISTA = 1
      do ISYM=1,nIrrep
        NA = NASH(ISYM)
        NO = NOSH(ISYM)
        if (NA /= 0) call SG_SSOTRA(iState,ISYM,LSM,NA,NO,TRA(ISTA),NCO,CI,TMP)
        ISTA = ISTA+NO**2
      end do
      call mma_deallocate(TMP)

    Case Default
      Write (u6,'(A)') ' RASSI_CITRA: Illegal WFTP value'
      Write (u6,'(2A)')'WFTP=',WFTP
      Call Abend()

  End Select
# ifdef _DEBUGPRINT_
    write(u6,*) ' DONE in  CITRA. norm=',ddot_(NCO,CI,1,CI,1)
# endif

  contains

  ! Performs the single-orbital transformations
  subroutine SG_SSOTRA(iState,ISYM,LSM,NA,NO,TRA,NCO,CI,TMP)

  use sguga, only: SGS, CIS, EXS
  use sguga, only: sg_epq_psi
  implicit none
  integer(kind=iwp), intent(in) :: iState, ISYM, LSM, NA, NO, NCO
  real(kind=wp), intent(in) :: TRA(NO,NO)
  real(kind=wp), intent(inout) :: CI(NCO)
  real(kind=wp), intent(out) :: TMP(NCO)
  integer(kind=iwp) :: IK, IKLEV, IL, IP, IPLEV, NI
  real(kind=wp) :: CKK, CPK, X
  integer(kind=iwp), allocatable :: ILEV(:)

  ! ILEV(IORB)=GUGA LEVEL CORRESPONDING TO A SPECIFIC ACTIVE ORBITAL
  ! OF SYMMETRY ISYM.

  call mma_allocate(ILEV,NA,Label='ILEV')

  NI = NO-NA
  IL = 0
  do IP=1,NA
    do
      IL = IL+1
      if (SGS(istate)%ISM(IL) == ISYM) exit
    end do
    ILEV(IP) = IL
  end do
  !Call RecPrt('SG_SSOTRA: CI',' ',CI,1,NCO)
  !Call RecPrt('TRA',' ',Tra,NO,NO)
  do IK=1,NA

  !   Write (6,*) 'ik=',ik

    IKLEV = ILEV(IK)
    TMP(:) = Zero
    do IP=1,NA
      IPLEV = ILEV(IP)
      CPK = TRA(NI+IP,NI+IK)
      if (IP == IK) CPK = CPK-One
      X = Half*CPK
      if (abs(X) < 1.0e-14_wp) cycle
      call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),IPLEV,IKLEV,X,LSM,CI,TMP)
    end do
  !  Call RecPrt('TMP',' ',TMP,1,NCO)

    CKK = TRA(NI+IK,NI+IK)
    X = Three-CKK
  !  Write (6,*) 'CKK, X=',CKK, X
    CI(:) = CI(:)+X*TMP(:)
  !  Call RecPrt('CI',' ',CI,1,NCO)

    do IP=1,NA
      IPLEV = ILEV(IP)
      CPK = TRA(NI+IP,NI+IK)
      if (IP == IK) CPK = CPK-One
      if (abs(CPK) < 1.0e-14_wp) cycle
      call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),IPLEV,IKLEV,CPK,LSM,TMP,CI)
    end do
  !  Call RecPrt('CI',' ',CI,1,NCO)

  end do
  call mma_deallocate(ILEV)

  end subroutine SG_SSOTRA

end subroutine SG_CITRA

subroutine SYG2SGU(IMODE,iState,LSYM,ICNFTAB,ISPNTAB,CIOLD,CINEW,IPGLOB)
  ! iState    : the split graph and CI array structure of SGUGA state iState
  ! IPGLOB    : print level
  ! IMODE=0 transforms a Symmetric Group CI array to SGUGA
  ! IMODE=1 transforms a Split GUGA CI array to Symm Group
  ! ...Configuration and Spin Coupling tables, fill this in later.
  ! CIOLD and CINEW are obvious.

  use sguga, only: CIS_All=>CIS, SGS_All=>SGS
  use Molcas, only: MxLev
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: One
  use Definitions, only: wp, iwp, u6


  implicit none
  integer(kind=iwp), intent(in) :: IMODE, iState, LSYM, ICNFTAB(*), ISPNTAB(*), IPGLOB
  real(kind=wp), intent(in) :: CIOLD(*)
  real(kind=wp), intent(_OUT_) :: CINEW(*)
  integer(kind=iwp), parameter :: IFUP2CS(0:1) = [2,1], MXCPI = 15, NBUFFER = 600
  integer(kind=iwp) :: I, ICASE(400), ICNF, ICNUM(NBUFFER), ICPL, ICSPLT, ICSYMG, IEL, IEL1, IEL2, IFORM, IFUP, IOCC, IORB, IREST, &
                       IWLKPOS, IWORD, IWRD, KCNF, KCNFINF, KCPL, KGSLIM, KGSORB, KSPNINF, KWALK(NBUFFER), LEV, MAXOP, MINOP, &
                       MIPWLK, MXWLK, NACTEL, NAPART, NCLSD, NCNF, NCONF, NCPL, NCSYMG, NHEAD, NLEV, NOCC, NODD, NOPEN, NSYM, &
                       NWALK, NWLKLST, NWRD
  real(kind=wp) :: PHASE(NBUFFER), PHS
  integer(kind=iwp), allocatable :: MWS2W(:), OrbArr(:)

  associate (SGS=>SGS_All(iState),CIS=>CIS_All(iState))

  ! Dereference CIS and SGS       for some data:
  NCONF = CIS%NCSF(LSYM)
  NWALK = CIS%nWalk
  call mma_allocate(MWS2W,NWALK,Label='MWS2W')
  NSYM = ICNFTAB(7)
  call SG_MSTOW(SGS,CIS,MWS2W,NSYM)
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
  ! Max nr of Split-Graph CSF's, stored as walks in the
  ! buffer KWALK
  MXWLK = NBUFFER/MIPWLK
  ! Pick up various data from conf and spin tables
  ! Unbutton Configuration table:

  MINOP = ICNFTAB(5)
  MAXOP = ICNFTAB(6)
  NSYM = ICNFTAB(7)
  if (LSYM /= ICNFTAB(8)) call abend()
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
  ! NCSYMG=Nr of Symmetric-Group CSF's treated so far.
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
        ORBARR(1:NOCC) = ICNFTAB(KCNF+NWRD*(ICNF-1):KCNF+NWRD*(ICNF-1)+NOCC-1)
        ! Loop over spin couplings
        do ICPL=1,NCPL
          ICASE(1:NLEV) = 0
          do I=1,NCLSD
            ICASE(ORBARR(I)) = 3
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
            NODD = (2*NODD-1)*I*(I-3)/2+NODD
            !if (NODD == 1) then
            !  IF(I == 2) PHS = -PHS
            !  IF(I == 3) PHS = -PHS
            !end if
            PHS = real((3+NODD*I*(I-1)*(2*I-7))/3,kind=wp)*PHS
          end do
          ! Pack the walk and add it to the list.
          call SG_PKWLK(NLEV,MIPWLK,1,KWALK(IWLKPOS),ICASE)
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
                ICSPLT = ICNUM(I)
                ICSYMG = NCSYMG+I
                CINEW(ICSPLT) = CIOLD(ICSYMG)*PHASE(I)
              end do
            else
              do I=1,NWLKLST
                ICSPLT = ICNUM(I)
                ICSYMG = NCSYMG+I
                CINEW(ICSYMG) = CIOLD(ICSPLT)*PHASE(I)
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
          ICASE(1:NLEV) = 0
          do I=1,NCLSD
            ICASE(ORBARR(I)) = 3
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
            NODD = (2*NODD-1)*I*(I-3)/2+NODD
            !if (NODD == 1) then
            !  if (I == 2) PHS = -PHS
            !  if (I == 3) PHS = -PHS
            !end if
            PHS = real((3+NODD*I*(I-1)*(2*I-7))/3,kind=wp)*PHS
          end do
          ! Pack the walk and add it to the list.
          call SG_PKWLK(NLEV,MIPWLK,1,KWALK(IWLKPOS),ICASE)
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
                ICSPLT = ICNUM(I)
                ICSYMG = NCSYMG+I
                CINEW(ICSPLT) = CIOLD(ICSYMG)*PHASE(I)
              end do
            else
              do I=1,NWLKLST
                ICSPLT = ICNUM(I)
                ICSYMG = NCSYMG+I
                CINEW(ICSYMG) = CIOLD(ICSPLT)*PHASE(I)
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
          ICASE(1:NLEV) = 0
          do I=1,NCLSD
            ICASE(ORBARR(I)) = 3
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
            NODD = (2*NODD-1)*I*(I-3)/2+NODD
            !if (NODD == 1) then
            !  if (I == 2) PHS = -PHS
            !  if (I == 3) PHS = -PHS
            !end if
            PHS = real((3+NODD*I*(I-1)*(2*I-7))/3,kind=wp)*PHS
          end do
          ! Pack the walk and add it to the list.
          call SG_PKWLK(NLEV,MIPWLK,1,KWALK(IWLKPOS),ICASE)
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
                ICSPLT = ICNUM(I)
                ICSYMG = NCSYMG+I
                CINEW(ICSPLT) = CIOLD(ICSYMG)*PHASE(I)
              end do
            else
              do I=1,NWLKLST
                ICSPLT = ICNUM(I)
                ICSYMG = NCSYMG+I
                CINEW(ICSYMG) = CIOLD(ICSPLT)*PHASE(I)
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
          ICASE(1:NLEV) = 0
          do I=1,NCLSD
            ICASE(ORBARR(I)) = 3
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
            NODD = (2*NODD-1)*I*(I-3)/2+NODD
            !if (NODD == 1) then
            !  if (I == 2) PHS = -PHS
            !  if (I == 3) PHS = -PHS
            !end if
            PHS = real((3+NODD*I*(I-1)*(2*I-7))/3,kind=wp)*PHS
          end do
          ! Pack the walk and add it to the list.
          call SG_PKWLK(NLEV,MIPWLK,1,KWALK(IWLKPOS),ICASE)
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
                ICSPLT = ICNUM(I)
                ICSYMG = NCSYMG+I
                CINEW(ICSPLT) = CIOLD(ICSYMG)*PHASE(I)
              end do
            else
              do I=1,NWLKLST
                ICSPLT = ICNUM(I)
                ICSYMG = NCSYMG+I
                CINEW(ICSYMG) = CIOLD(ICSPLT)*PHASE(I)
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
      ICSPLT = ICNUM(I)
      ICSYMG = NCSYMG+I
      CINEW(ICSPLT) = CIOLD(ICSYMG)*PHASE(I)
    end do
  else
    do I=1,NWLKLST
      ICSPLT = ICNUM(I)
      ICSYMG = NCSYMG+I
      CINEW(ICSYMG) = CIOLD(ICSPLT)*PHASE(I)
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

  end associate

end subroutine SYG2SGU

subroutine SG_MSTOW(SGS,CIS,MWS2W,nSym)
  ! Purpose: From the list of packed up- and downwalks, construct
  ! the table MWS2W, such that MAW sums can be translated to the
  ! corresponding walks of the Split-GUGA.

  use sguga, only: CIStruct, SGStruct
  use stdalloc, only: mma_allocate, mma_deallocate
  use Definitions, only: iwp


  implicit none
  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
  integer(kind=iwp), intent(_OUT_) :: MWS2W(*)
  integer(kind=iwp), intent(in) :: nSym
  integer(kind=iwp) :: IC, IDOFF, IDV, IDW, IDWTOT, ISYDWN, ISYUP, IUOFF, IUV, IUW, IUWTOT, LEV, MS, MV, NDWN, NUP
  integer(kind=iwp), allocatable :: ICS(:)

  call mma_allocate(ICS,SGS%nLev,Label='ICS')

  do MV=1,CIS%nMidV
    do ISYUP=1,NSYM
      NUP = CIS%NOW(1,ISYUP,MV)
      if (NUP == 0) cycle
      IUOFF = CIS%IOW(1,ISYUP,MV)/CIS%nIpWlk
      do IUW=1,NUP
        IUWTOT = IUOFF+IUW
        ! Unpack upper walk to ICS()
        call SG_UPKWLK(SGS%nLev-SGS%MidLev,CIS%nIpWlk,1,CIS%ICase(1+CIS%nIpWlk*(IUWTOT-1)),ICS(SGS%MidLev+1))
        MS = 0
        IUV = 1
        do LEV=SGS%nLev,SGS%MidLev+1,-1
          IC = ICS(LEV)
          MS = MS+SGS%MAW(IUV,IC)
          IUV = SGS%DOWN(IUV,IC)
        end do
        MWS2W(MS) = IUWTOT
      end do
    end do
  end do

  do MV=1,CIS%nMidV
    do ISYDWN=1,NSYM
      NDWN = CIS%NOW(2,ISYDWN,MV)
      if (NDWN == 0) cycle
      IDOFF = CIS%IOW(2,ISYDWN,MV)/CIS%nIpWlk
      do IDW=1,NDWN
        IDWTOT = IDOFF+IDW
        ! Unpack lower walk to ICS()
        call SG_UPKWLK(SGS%MidLev,CIS%nIpWlk,1,CIS%ICase(1+CIS%nIpWlk*(IDWTOT-1)),ICS)
        MS = 0
        IDV = SGS%nVert
        do LEV=1,SGS%MidLev
          IC = ICS(LEV)
          IUV = SGS%UP(IDV,IC)
          MS = MS+SGS%MAW(IUV,IC)
          IDV = IUV
        end do
        MWS2W(MS) = IDWTOT
      end do
    end do
  end do

  call mma_deallocate(ICS)

end subroutine SG_MSTOW

subroutine W2SGORD(SGS,CIS,MWS2W,NLIST,KWALK,ICNUM)
  ! Purpose: Given a list of bit-packed total walks,
  ! translate this into a list of elements of a CI array.
  ! MXCPI= Max nr of case numbers packed in one integer.
  !
  ! For a wave function in Split GUGA storage structure,
  ! given KWALK(J,I) with J=1..MIPWLK that contains the
  ! complete Guga walk, as a packed array of case numbers, construct
  ! ICNUM(I), which is the index of this configuration in the
  ! Split-Guga scheme.
  ! Use MAWS to WLK table, MWS2W

  use sguga, only: CIStruct, SGStruct
  use Symmetry_Info, only: MUL
  use stdalloc, only: mma_allocate, mma_deallocate
  use Definitions, only: iwp

  implicit none
  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
  integer(kind=iwp), intent(in) :: MWS2W(*), NLIST, KWALK(*)
  integer(kind=iwp), intent(out) :: ICNUM(NLIST)
  integer(kind=iwp) :: IC, ICONF, IDV, IDW, IOFF, ISYCI, ISYDWN, ISYUP, IUV, IUW, LDIM, LEV, MAWSD, MAWSU, MIPWLK, MV
  integer(kind=iwp), allocatable :: ICS(:)
  integer(kind=iwp), parameter :: MXCPI = 15

  ! Nr of integers used to store each total walk:
  MIPWLK = 1+(SGS%nLev-1)/MXCPI
  ! Allocate scratch space for case numbers:
  call mma_allocate(ICS,SGS%nLev,Label='ICS')

  do ICONF=1,NLIST
    ! Unpack total walk to ICS()
    call SG_UPKWLK(SGS%nLev,MIPWLK,1,KWALK((ICONF-1)*MIPWLK+1),ICS)

    ! Follow upper walk down to SGS%MidLev:
    MAWSU = 0
    IUV = 1
    ISYUP = 1
    IDV = -1000000000
    do LEV=SGS%nLev,SGS%MidLev+1,-1
      IC = ICS(LEV)
      if (((IC+1)/2) == 1) ISYUP = MUL(SGS%ISM(LEV),ISYUP)
      IDV = SGS%DOWN(IUV,IC)
      MAWSU = MAWSU+SGS%MAW(IUV,IC)
      IUV = IDV
    end do
    ! We have found the midvertex number:
    MV = IDV+1-SGS%MVSta
    ! Follow lower walk down to SGS%MidLev:
    ISYDWN = 1
    MAWSD = 0
    do LEV=SGS%MidLev,1,-1
      IC = ICS(LEV)
      if (((IC+1)/2) == 1) ISYDWN = MUL(SGS%ISM(LEV),ISYDWN)
      IDV = SGS%DOWN(IUV,IC)
      MAWSD = MAWSD+SGS%MAW(IUV,IC)
      IUV = IDV
    end do
    IUW = MWS2W(MAWSU)
    IDW = MWS2W(MAWSD)
    ! Subtract the offsets:
    IUW = IUW-CIS%IOW(1,ISYUP,MV)/CIS%nIpWlk
    IDW = IDW-CIS%IOW(2,ISYDWN,MV)/CIS%nIpWlk
    ! Split-Guga storage scheme: an element in a set of matrices.
    ! Offset to the matrix we want is CIS%IOCSF(ISYUP,MV,ISYCI).
    ! Leading dimension=nr of upwalks in this block.
    ISYCI = MUL(ISYUP,ISYDWN)
    IOFF = CIS%IOCSF(ISYUP,MV,ISYCI)
    LDIM = CIS%NOW(1,ISYUP,MV)
    ICNUM(ICONF) = IOFF+IUW+LDIM*(IDW-1)
  end do

  call mma_deallocate(ICS)

end subroutine W2SGORD

subroutine SG_PKWLK(N,IPWLK,NWALK,IWALK,ICASE)
  ! PURPOSE: PACK THE GUGA STEP NUMBERS INTO THE ARRAY IWALK.
  ! EACH OF THE NWALK WALKS HAS N STEP NUMBERS, 2 BITS EACH,
  ! AT MOST nPack TO AN INTEGER ELEMENT OF IWALK, EACH NEW WALK
  ! IS ALIGNED ON INTEGERS.
  ! NOTE: Can be used for upper, lower, or total walks, so the
  ! number of integers used for each walk (IPWLK) is given as
  ! call parameter.

  use sguga, only: nPack
  use Definitions, only: iwp


  implicit none
  integer(kind=iwp), intent(in) :: N, IPWLK, NWALK, ICASE(N,NWALK)
  integer(kind=iwp), intent(_OUT_) :: IWALK(*)
  integer(kind=iwp) :: I, IPOS, IWORD, J, L, LEND, LSTA

  IPOS = 0
  do I=1,NWALK
    LEND = 0
    do J=1,IPWLK
      LSTA = LEND+1
      LEND = min(LSTA+nPack-1,N)
      IPOS = IPOS+1
      IWORD = 0
      do L=LEND,LSTA,-1
        IWORD = 4*IWORD+ICASE(L,I)
      end do
      IWALK(IPOS) = IWORD
    end do
  end do

end subroutine SG_PKWLK

subroutine SG_UPKWLK(N,IPWLK,NWALK,IWALK,ICASE)

  use sguga, only: nPack
  use Definitions, only: iwp

  implicit none
  integer(kind=iwp), intent(in) :: N, IPWLK, NWALK, IWALK(*)
  integer(kind=iwp), intent(out) :: ICASE(N,NWALK)
  integer(kind=iwp) :: I, IPOS, IWORD, J, L, LEND, LSTA, NEXT

  ! See companion subroutine PKWLK.
  IPOS = 0
  do I=1,NWALK
    LEND = 0
    do J=1,IPWLK
      LSTA = LEND+1
      LEND = min(LSTA+nPack-1,N)
      IPOS = IPOS+1
      IWORD = IWALK(IPOS)
      do L=LSTA,LEND
        NEXT = IWORD/4
        ICASE(L,I) = IWORD-4*NEXT
        IWORD = NEXT
      end do
    end do
  end do

end subroutine SG_UPKWLK

end module SGUGA_Transform
