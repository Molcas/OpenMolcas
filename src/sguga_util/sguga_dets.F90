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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!               2026, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
!#define _TEST_CSFDET_
! CSF <-> determinant expansion of the open shells (CDC, Mk_CDC) and the
! CSF/determinant reordering tables (Mk_SDREO_SGUGA).
! (Default public, as SGUGA_Core; the interface is the module SGUGA.)
module SGUGA_Dets

use SGUGA_Core

#include "intent.fh"

implicit none

contains

subroutine Mk_CDC(SGS)

  type(SGStruct), intent(inout) :: SGS

  integer(kind=iwp) :: iOpen, mult, ms2, nDet, nCSF, MinOpen, MaxOpen, iSym, iDOC
  integer(kind=iwp) :: nClosed, nClosedOCP, nSOC, iSOC, nDOC, iOC
  integer(kind=iwp) :: NoNew, nConf, nElec, nDetSC
  integer :: i, j
  integer :: iStart
  integer(kind=iwp), allocatable :: nSOCSym(:)
  integer(kind=iwp), allocatable :: detspin(:,:)
  integer(kind=iwp), allocatable :: twos(:,:)

  integer :: DOCOcc(SGS%nLev)
  integer :: SOCOcc(SGS%nLev)
  integer(kind=iwp) :: nDOCList
  integer(kind=iwp) :: nSOCList

  allocate(SGS%CDCInfo(0:SGS%nLev,SGS%iSpin:SGS%iSpin))
  SGS%CDCInfo(:,:)%nDet = 0

  mult = SGS%iSpin

  MinOpen = mult-1
  MaxOpen = min(SGS%nActEl,2*SGS%nLev-SGS%nActel)
  if (mod(MaxOpen-MinOpen,2) /= 0) MaxOpen = MaxOpen-1

  call mma_allocate(nSOCSym,SGS%nSym,Label='nSOCSym')
  do iOpen=MinOpen,MaxOpen,2

    ms2 = mult-1

    call count_detspin(iOpen,ms2,nDet)
    call count_twos(iOpen,mult,nCSF)

    SGS%CDCInfo(iOpen,mult)%nDet = nDet
    SGS%CDCInfo(iOpen,mult)%nCSF = nCSF

    nClosed = (SGS%nActEl-iOpen)/2
    SGS%CDCInfo(iOpen,Mult)%nClosed = nClosed

    nClosedOCP = binomial(nClosed,SGS%nLev)
    SGS%CDCInfo(iOpen,Mult)%nClosedOCP = nClosedOCP

    nSOC = binomial(iOpen,SGS%nLev)

    if (nClosed > 0) then
      call mma_allocate(SGS%CDCInfo(iOpen,mult)%DOC, &
                        nClosed, &
                        nClosedOCP, &
                        Label='DOC')
      iDOC = 0
      call Make_DOCs(1,1)
      if (SG_Checks) then  ! sanity check, only when verifying CI developments
        if (iDOC /= SGS%CDCInfo(iOpen,mult)%nClosedOCP) then
          write(u6,*) 'DOC generation failed'
          write(u6,*) 'Generated = ',iDOC
          write(u6,*) 'Expected  = ',SGS%CDCInfo(iOpen,mult)%nClosedOCP
          call Abend()
        end if
      end if
    end if

    nSOCSym(:) = 0

    call mma_allocate(SGS%CDCInfo(iOpen,mult)%nOCP, &
                      SGS%nSym, &
                      Label='nOCP')

    call mma_allocate(SGS%CDCInfo(iOpen,mult)%nOpenOCP, &
                      SGS%nSym, &
                      Label='nOpenOCP')
    SGS%CDCInfo(iOpen,mult)%nOpenOCP(:) = 0

    call Count_OpenShell_OCPs(SGS%ISm, &
                              SGS%nLev, &
                              iOpen, &
                              SGS%nSym, &
                              SGS%CDCInfo(iOpen,mult)%nOpenOCP)

    SGS%CDCInfo(iOpen,mult)%nOCP(:) = SGS%CDCInfo(iOpen,mult)%nOpenOCP(:)*SGS%CDCInfo(iOpen,mult)%nClosedOCP

    if (iOpen > 0) then
      call mma_allocate(SGS%CDCInfo(iOpen,mult)%SOC, &
                        iOpen, &
                        nSOC, &
                        Label='SOC')

      call mma_allocate(SGS%CDCInfo(iOpen,mult)%SOCSym, &
                        nSOC, &
                        Label='SOCSym')
      iSOC = 0
      nSOCSym(:) = 0
      call Make_SOCs(1,1,1)

      if (SG_Checks) then  ! sanity check, only when verifying CI developments
        if (iSOC /= nSOC) then
          write(u6,*) 'SOC generation failed'
          write(u6,*) 'Generated = ',iSOC
          write(u6,*) 'Expected  = ',nSOC
          call Abend()
        end if
        do iSym=1,SGS%nSym
          if (nSOCSym(iSym) /= SGS%CDCInfo(iOpen,mult)%nOpenOCP(iSym)) then
            write(u6,*) 'SOC symmetry count mismatch'
            write(u6,*) 'ISYM     = ',iSym
            write(u6,*) 'GENERATED= ',nSOCSym(iSym)
            write(u6,*) 'EXPECTED = ',SGS%CDCInfo(iOpen,mult)%nOpenOCP(iSym)
            call Abend()
          end if
        end do
      end if

    end if

    ! Allocate and initialize the reverse list: iDOC,iSOC -> iOC

    if (allocated(SGS%CDCInfo(iOpen,mult)%DOC)) then
      nDOCList = size(SGS%CDCInfo(iOpen,mult)%DOC,2)
    else
      nDOCList = 1
    end if
    if (allocated(SGS%CDCInfo(iOpen,mult)%SOC)) then
      nSOCList = size(SGS%CDCInfo(iOpen,mult)%SOC,2)
    else
      nSOCList = 1
    end if
    call mma_allocate(SGS%CDCInfo(iOpen,mult)%OCReverse, &
                      nDOCList,nSOCList,SGS%nSym, &
                      Label='OCReverse')

    SGS%CDCInfo(iOpen,mult)%OCReverse = 0_iwp

    if (SGS%CDCInfo(iOpen,mult)%nDet > 0 .and. &
        SGS%CDCInfo(iOpen,mult)%nCSF > 0) then

      call generate_detspin(iOpen,mult-1,detspin,nDet)

      call mma_allocate( &
        SGS%CDCInfo(iOpen,mult)%DetSpin, &
        [lbound(detspin,1),ubound(detspin,1)], &
        [lbound(detspin,2),ubound(detspin,2)], &
        Label='DetSpin')

      SGS%CDCInfo(iOpen,mult)%DetSpin(:,:) = detspin(:,:)

      call generate_twos(iOpen,mult,twos,nCSF)
      call mma_allocate( &
        SGS%CDCInfo(iOpen,mult)%Twos, &
        [lbound(twos,1),ubound(twos,1)], &
        [lbound(twos,2),ubound(twos,2)], &
        Label='Twos')

      SGS%CDCInfo(iOpen,mult)%Twos(:,:) = twos(:,:)

      call mma_allocate( &
        SGS%CDCInfo(iOpen,mult)%CDC, &
        nDet,nCSF, &
        Label='CDC')

      call sg_csfdet( &
        iOpen, &
        nDet, &
        detspin, &
        nCSF, &
        twos, &
        SGS%CDCInfo(iOpen,mult)%CDC)

! keep this for the time being.
      if (mult == 1) then

        nDetSC = max(1,nDet/2)
        SGS%CDCInfo(iOpen,mult)%nDetSC = nDetSC

        call mma_allocate( &
          SGS%CDCInfo(iOpen,mult)%DTOC, &
          nDetSC, &
          nCSF, &
          Label='DTOC')

        call mma_allocate( &
          SGS%CDCInfo(iOpen,mult)%Pair, &
          nDet, &
          Label='Pair')

        if (nDet > 1) then
          call mk_dtoc_singlet( &
            nDet, &
            nCSF, &
            SGS%CDCInfo(iOpen,mult)%DetSpin, &
            SGS%CDCInfo(iOpen,mult)%CDC, &
            SGS%CDCInfo(iOpen,mult)%Pair, &
            SGS%CDCInfo(iOpen,mult)%DTOC)
        else if (nDet == 1) then
          SGS%CDCInfo(iOpen,mult)%Pair(1) = 1
          SGS%CDCInfo(iOpen,mult)%DTOC(1,1) = SGS%CDCInfo(iOpen,mult)%CDC(1,1)
        else
          write(u6,*) 'Mk_CDC: illegal nDet value: nDet<0'
          call Abend()
        end if
      end if

      call mma_deallocate(detspin)
      call mma_deallocate(twos)

    end if

    if (SG_Checks) then  ! sanity check, only when verifying CI developments
      do iSym=1,SGS%nSym
        if (SGS%CDCInfo(iOpen,mult)%nOCP(iSym) /= &
            SGS%CDCInfo(iOpen,mult)%nOpenOCP(iSym)* &
            SGS%CDCInfo(iOpen,mult)%nClosedOCP) then
          call Abend()
        end if
      end do

      if (sum(SGS%CDCInfo(iOpen,mult)%nOpenOCP) /= binomial(iOpen,SGS%nLev)) then
        call Abend()
      end if
    end if

    ! Compute the number of orbital occupatation configurations per symmetry.
    ! Here the order doesn't matter, we will deal with that in the next loop
    ! where we copy the Lucia order.

    call mma_allocate(SGS%CDCInfo(iOpen,Mult)%nOCSym,SGS%nSym,Label='nOCSym')

    SGS%CDCInfo(iOpen,mult)%nOCSym(:) = &
      SGS%CDCInfo(iOpen,mult)%nClosedOCP* &
      SGS%CDCInfo(iOpen,mult)%nOpenOCP(:)

    call mma_allocate(SGS%CDCInfo(iOpen,mult)%OCDOC, &
                      maxval(SGS%CDCInfo(iOpen,mult)%nOCSym), &
                      SGS%nSym, &
                      Label='OCDOC')

    call mma_allocate(SGS%CDCInfo(iOpen,mult)%OCSOC, &
                      maxval(SGS%CDCInfo(iOpen,mult)%nOCSym), &
                      SGS%nSym, &
                      Label='OCSOC')

    SGS%CDCInfo(iOpen,mult)%OCDOC(:,:) = 0
    SGS%CDCInfo(iOpen,mult)%OCSOC(:,:) = 0
  end do

  call mma_deallocate(nSOCSym)

  call mma_allocate(SGS%OC,SGS%nActel,Label='OC')

! Now loop over the Orbital Configurations a la Lucia.
! For each OC identify the symmtry and store away the iDOC and iSOC indicies.

! Reset and use now as a counter
  do iOpen=MinOpen,MaxOpen,2
    SGS%CDCInfo(iOpen,mult)%nOCSym(:) = 0
  end do

  nElec = SGS%nActEl
  nConf = 0
  NoNew = 0
  call Nxt_OC(SGS%OC,nElec,SGS%nLev,1,NoNew)

  do while (NoNew == 0)
    nConf = nConf+1

    call Decode_OC(SGS%OC,nElec, &
                   nDOC,DOCOcc, &
                   nSOC,SOCOcc)

    if (nSOC < MinOpen .or. nSOC > MaxOpen) then
      call Nxt_OC(SGS%OC,nElec,SGS%nLev,0,NoNew)
      cycle
    end if

    iDOC = 1
    do i=1,nDOC
      if (i == 1) then
        iStart = 1
      else
        iStart = DOCOcc(i-1)+1
      end if
      do j=iStart,DOCOcc(i)-1
        iDOC = iDOC+Binomial(nDOC-i,SGS%nLev-j)
      end do
    end do

    if (SG_Checks) then  ! sanity check, only when verifying CI developments
      if (nDOC > 0) then
        do i=1,nDOC
          if (DOCOcc(i) /= SGS%CDCInfo(nSOC,mult)%DOC(i,iDOC)) then
            write(u6,*) 'DOC ranking mismatch'
            write(u6,'(A,20I4)') 'DOC OCC = ',DOCOcc(1:nDOC)
            write(u6,'(A,I6)') 'iDOC    = ',iDOC
            write(u6,'(A,20I4)') 'DOC TAB = ', &
              SGS%CDCInfo(nSOC,mult)%DOC(1:nDOC,iDOC)
            call Abend()
          end if
        end do
      end if
    end if

    iSOC = 1
    do i=1,nSOC

      if (i == 1) then
        iStart = 1
      else
        iStart = SOCOcc(i-1)+1
      end if

      do j=iStart,SOCOcc(i)-1
        iSOC = iSOC+Binomial(nSOC-i,SGS%nLev-j)
      end do

    end do

    if (SG_Checks) then  ! sanity check, only when verifying CI developments
      if (nSOC > 0) then
        do i=1,nSOC
          if (SOCOcc(i) /= SGS%CDCInfo(nSOC,mult)%SOC(i,iSOC)) then
            write(u6,*) 'SOC ranking mismatch'
            write(u6,'(A,20I4)') 'SOC OCC = ',SOCOcc(1:nSOC)
            write(u6,'(A,I6)') 'iSOC    = ',iSOC
            write(u6,'(A,20I4)') 'SOC TAB = ', &
              SGS%CDCInfo(nSOC,mult)%SOC(1:nSOC,iSOC)
            call Abend()
          end if
        end do
      end if
!  Write(u6,'(A,4I8)') 'OC:', nConf,nSOC,iDOC,iSOC
    end if

    if (nSOC == 0) then
      iSym = 1
      iSOC = 1
    else
      iSym = SOC_Sym(SGS%ISm,nSOC,SOCOcc)

      if (SG_Checks) then  ! sanity check, only when verifying CI developments
        do i=1,nSOC
          if (SOCOcc(i) /= SGS%CDCInfo(nSOC,mult)%SOC(i,iSOC)) then
            write(u6,*) 'SOC ranking mismatch'
            write(u6,'(A,20I4)') 'SOC OCC = ',SOCOcc(1:nSOC)
            write(u6,'(A,I6)') 'iSOC    = ',iSOC
            write(u6,'(A,20I4)') 'SOC TAB = ', &
              SGS%CDCInfo(nSOC,mult)%SOC(1:nSOC,iSOC)
            call Abend()
          end if
        end do

        if (iSym /= SGS%CDCInfo(nSOC,mult)%SOCSym(iSOC)) then
          write(u6,*) 'SOC symmetry mismatch'
          call Abend()
        end if
      end if
    end if

    SGS%CDCInfo(nSOC,mult)%nOCSym(iSym) = SGS%CDCInfo(nSOC,mult)%nOCSym(iSym)+1

    iOC = SGS%CDCInfo(nSOC,mult)%nOCSym(iSym)

    SGS%CDCInfo(nSOC,mult)%OCDOC(iOC,iSym) = iDOC
    SGS%CDCInfo(nSOC,mult)%OCSOC(iOC,iSym) = iSOC

    SGS%CDCInfo(nSOC,mult)%OCReverse(iDOC,iSOC,iSym) = iOC

    call Nxt_OC(SGS%OC,nElec,SGS%nLev,0,NoNew)

  end do

contains

  subroutine mk_dtoc_singlet(nDet,nCSF,DetSpin,CDC,Pair,DTOC)

    use Definitions, only: wp, iwp
    implicit none

    integer(iwp), intent(in) :: nDet
    integer(iwp), intent(in) :: nCSF

    integer(iwp), intent(in) :: DetSpin(:,:)

    real(wp), intent(in) :: CDC(nDet,nCSF)

    integer(iwp), intent(out) :: Pair(nDet)

    real(wp), intent(out) :: DTOC(nDet/2,nCSF)

    integer(iwp) :: iDet
    integer(iwp) :: jDet
    integer(iwp) :: iRow

    Pair(:) = 0

!
! Determine spin-reversal partner for every determinant
!
    do iDet=1,nDet

      do jDet=1,nDet

        if (all(DetSpin(:,iDet) == 1-DetSpin(:,jDet))) then
          Pair(iDet) = jDet
          exit
        end if

      end do

      if (Pair(iDet) == 0) then
        write(6,*) 'MK_DTOC_SINGLET: partner not found'
        write(6,*) 'IDET =',iDet
        call Abend()
      end if

    end do

!
! Generate Lucia DTOC matrix.
!
! Keep one determinant from every spin-reversal pair
! and renormalize by sqrt(2).
!
    iRow = 0

    do iDet=1,nDet

      ! Lucia/SPNCOM ordering:
      ! keep determinants with first spin = alpha

      if (DetSpin(1,iDet) == 1) then

        iRow = iRow+1

        DTOC(iRow,:) = sqrt(2.0_wp)*CDC(iDet,:)

      end if

    end do

    if (iRow /= max(1,nDet/2)) then
      write(6,*) 'MK_DTOC_SINGLET: incorrect number of pairs'
      write(6,*) 'Generated =',iRow
      write(6,*) 'Expected  =',max(1,nDet/2)
      call Abend()
    end if

  end subroutine mk_dtoc_singlet

  recursive subroutine Make_SOCs(iPos,iStart,iSym)

    integer(iwp), intent(in) :: iPos
    integer(iwp), intent(in) :: iStart
    integer(iwp), intent(in) :: iSym

    integer(iwp) :: i
    integer(iwp), save :: work(MxLev)

    if (iPos > iOpen) then

      iSOC = iSOC+1

      SGS%CDCInfo(iOpen,mult)%SOC(:,iSOC) = work(1:iOpen)
      SGS%CDCInfo(iOpen,mult)%SOCSym(iSOC) = iSym

      nSOCSym(iSym) = nSOCSym(iSym)+1

#ifdef _DEBUGPRINT_
      write(u6,*) 'ISm = ',SGS%ISm(1:SGS%nLev)
      write(u6,*) 'NOPEN=',iOpen, &
        ' SOC=',iSOC, &
        ' OCC=',work(1:iOpen), &
        ' SYM=',iSym
#endif
      return
    end if

    do i=iStart,SGS%nLev-(iOpen-iPos)

      work(iPos) = i

      call Make_SOCs(iPos+1, &
                     i+1, &
                     Mul(iSym,SGS%ISm(i)))

    end do

  end subroutine Make_SOCs
  recursive subroutine Make_DOCs(iPos,iStart)

    integer(iwp), intent(in) :: iPos
    integer(iwp), intent(in) :: iStart

    integer(iwp) :: i
    integer(iwp), save :: work(MxLev)
!
!  Completed DOC
!
    if (iPos > nClosed) then
      iDOC = iDOC+1
      SGS%CDCInfo(iOpen,mult)%DOC(:,iDOC) = work(1:nClosed)
#ifdef _DEBUGPRINT_
      write(u6,'(A,I4,A,100I4)') 'NOPEN=',iOpen, &
        ' DOC=',iDOC, &
        SGS%CDCInfo(iOpen,mult)%DOC(:,iDOC)
#endif
      return
    end if
!
!  Lexical combinations
!
    do i=iStart,SGS%nLev-(nClosed-iPos)
      work(iPos) = i
      call Make_DOCs(iPos+1,i+1)
    end do

  end subroutine Make_DOCs

  subroutine count_detspin(nopen,ms2,ndet)

    use Definitions, only: iwp

    implicit none

    integer(iwp), intent(in) :: nopen
    integer(iwp), intent(in) :: ms2
    integer(iwp), intent(out) :: ndet

    integer(iwp) :: ibit
    integer(iwp) :: mx

    ndet = 0
    mx = 2**nopen

    do ibit=0,mx-1
      if (2*popcnt(ibit)-nopen == ms2) ndet = ndet+1
    end do

  end subroutine count_detspin

  subroutine count_twos(nopen,mult,ncsf)

    use Definitions, only: iwp

    implicit none

    integer(iwp), intent(in) :: nopen
    integer(iwp), intent(in) :: mult
    integer(iwp), intent(out) :: ncsf

    integer(iwp) :: ibit
    integer(iwp) :: i
    integer(iwp) :: ms2l
    integer(iwp) :: mx

    logical :: upper

    ncsf = 0
    mx = 2**nopen

    do ibit=0,mx-1

      ms2l = 0
      upper = .true.

      do i=1,nopen

        if (ibits(ibit,i-1,1) == 1) then
          ms2l = ms2l+1
        else
          ms2l = ms2l-1
        end if

        if (ms2l < 0) then
          upper = .false.
          exit
        end if

      end do

      if (.not. upper) cycle
      if (ms2l /= mult-1) cycle

      ncsf = ncsf+1

    end do

  end subroutine count_twos

end subroutine Mk_CDC

!***********************************************************************
! Construct the SD reordering information required by the SGUGA-based
! CAS SD-CSF codec.
!
! Initial version:
!   Loop over the alpha and beta string symmetry pairs compatible with
!   the requested total spatial symmetry.
!***********************************************************************

subroutine Mk_SDREO_SGUGA(SGS)

  use Definitions, only: iwp
  use Symmetry_Info, only: Mul

  implicit none

  type(SGStruct), intent(inout) :: SGS

  integer(kind=iwp) :: iRefSm
  integer(kind=iwp) :: iASm
  integer(kind=iwp) :: iBSm
  integer(kind=iwp) :: nSym
  integer(kind=iwp) :: iDet

  integer(kind=iwp) :: nAlpha
  integer(kind=iwp) :: nBeta
  integer(kind=iwp) :: nAlphaS
  integer(kind=iwp) :: nBetaS
  integer(kind=iwp) :: iAlphaS
  integer(kind=iwp) :: iBetaS

  integer(kind=iwp) :: nAlphaSSym(SGS%nSym)
  integer(kind=iwp) :: nBetaSSym(SGS%nSym)
  integer(kind=iwp), allocatable :: AlphaStrings(:,:)
  integer(kind=iwp), allocatable :: BetaStrings(:,:)
  integer(kind=iwp) :: iDetOC(SGS%nActEl)
  integer(kind=iwp) :: iDetSP(SGS%nActEl)
  integer(kind=iwp) :: iSign

  integer(kind=iwp) :: nOpen, i, iOpen, mult, iDOC, iSOC
  integer(kind=iwp) :: nDouble
  integer(kind=iwp) :: DOCOcc(SGS%nActel)
  integer(kind=iwp) :: SOCOcc(SGS%nActel)
  integer(kind=iwp) :: OpenSP(SGS%nActel)

  integer(kind=iwp) :: iOC
  integer(kind=iwp) :: nOC
  integer(kind=iwp) :: iSpinPattern
  integer(kind=iwp) :: iOpenOffset(0:SGS%nLev)
  integer(kind=iwp) :: jOpen
  integer(kind=iwp) :: MinOpen
  integer(kind=iwp) :: MaxOpen
  integer(kind=iwp) :: nOCSector
  integer(kind=iwp) :: nDetSector
  integer(kind=iwp) :: iAdr, nDetTotal
  integer(kind=iwp) :: iAdrSDConfOrder

  nSym = SGS%nSym
  mult = SGS%iSpin

  !
  ! There is one SD reorder list for each total symmetry.
  !
  allocate(SGS%SDREO(1:nSym))

  ! Highest-Ms component:
  !
  !   nAlpha + nBeta = nActEl
  !   nAlpha - nBeta = 2*Ms = iSpin - 1
  !
  nAlpha = (SGS%nActEl+SGS%iSpin-1_iwp)/2_iwp
  nBeta = (SGS%nActEl-SGS%iSpin+1_iwp)/2_iwp

  call Count_Strings_By_Symmetry(SGS,nAlpha,nAlphaSSym)
  call Count_Strings_By_Symmetry(SGS,nBeta,nBetaSSym)

  ! Outer loop over possible total spatial symmetries.
  do iRefSm=1,nSym

    MinOpen = mult-1_iwp

    MaxOpen = min(SGS%nActEl,2_iwp*SGS%nLev-SGS%nActEl)
    if (mod(MaxOpen-MinOpen,2_iwp) /= 0_iwp) MaxOpen = MaxOpen-1_iwp

    iOpenOffset(:) = 0_iwp
    iAdr = 0_iwp
    do jOpen=MinOpen,MaxOpen,2_iwp
      iOpenOffset(jOpen) = iAdr
      nOCSector = SGS%CDCInfo(jOpen,mult)%nOCSym(iRefSm)
      nDetSector = SGS%CDCInfo(jOpen,mult)%nDet

      iAdr = iAdr+nOCSector*nDetSector
    end do
    nDetTotal = iAdr

    if (allocated(SGS%SDREO(iRefSm)%List)) deallocate(SGS%SDREO(iRefSm)%List)

    call mma_allocate(SGS%SDREO(iRefSm)%List,nDetTotal,Label='List')

    SGS%SDREO(iRefSm)%List(:) = 0_iwp

    ! Loop over symmetry blocks in the same order as LUCIA.
    iDet = 0_iwp

    do iASm=1,nSym

      iBSm = Mul(iASm,iRefSm)

      nAlphaS = nAlphaSSym(iASm)
      nBetaS = nBetaSSym(iBSm)

      if (nAlphaS*nBetaS == 0) cycle

      call mma_allocate(AlphaStrings,max(1_iwp,nAlpha),nAlphaS,Label='AStrings')
      call mma_allocate(BetaStrings,max(1_iwp,nBeta),nBetaS,Label='BStrings')

      AlphaStrings(:,:) = 0_iwp
      BetaStrings(:,:) = 0_iwp

      call GenerateStringSet(nAlpha,iASm,nAlphaS,AlphaStrings)
      call GenerateStringSet(nBeta,iBSm,nBetaS,BetaStrings)

      do iBetaS=1,nBetaS
        do iAlphaS=1,nAlphaS

          call BuildOrderedDeterminant(AlphaStrings(:,iAlphaS),BetaStrings(:,iBetaS), &
                                       nAlpha,nBeta,iDetOC,iDetSP,iSign)

          call ExtractOccupationPattern(iDetOC,iDetSP,SGS%nActEl, &
                                        DOCOcc,SOCOcc,OpenSP,nDouble,nOpen)

          iOpen = nOpen

          if (SG_Checks) then  ! sanity check, only when verifying CI developments
            if (2_iwp*nDouble+iOpen /= SGS%nActEl) then
              write(u6,*) 'Mk_SDREO_SGUGA_new: inconsistent occupation counts'
              call Abend()
            end if
            if (mod(iOpen,2_iwp) /= mod(SGS%nActEl,2_iwp)) then
              write(u6,*) 'Mk_SDREO_SGUGA_new: invalid number of open orbitals'
              call Abend()
            end if

            if (iOpen > 0_iwp) then
              if (count(OpenSP(1:iOpen) == 1_iwp) /= nAlpha-nDouble) then
                write(u6,*) 'Mk_SDREO_SGUGA_new: incorrect open alpha count'
                call Abend()
              end if
              if (count(OpenSP(1:iOpen) == -1_iwp) /= nBeta-nDouble) then
                write(u6,*) 'Mk_SDREO_SGUGA_new: incorrect open beta count'
                call Abend()
              end if
            else
              if (nAlpha /= nDouble .or. nBeta /= nDouble) then
                write(u6,*) 'Mk_SDREO_SGUGA_new: invalid closed-shell determinant'
                call Abend()
              end if
            end if
            if (nDouble /= 0_iwp .and. .not. allocated(SGS%CDCInfo(nOpen,mult)%DOC)) then
              write(u6,*) 'Mk_SDREO_SGUGA_new: ...%DOC not allocated while nDouble /= 0!'
              call Abend()
            end if
          end if
          if (nDouble == 0_iwp) then
            iDOC = 1_iwp
          else
            call FindDOCIndex(iDetOC,nAlpha+nBeta, &
                              SGS%CDCInfo(nOpen,mult)%DOC,iDOC)
          end if

          if (nOpen == 0_iwp) then
            iSOC = 1_iwp
          else
            call FindSOCIndex(iDetOC,nAlpha+nBeta, &
                              SGS%CDCInfo(nOpen,mult)%SOC,iSOC)
          end if

          nOC = SGS%CDCInfo(nOpen,mult)%nOCSym(iRefSm)

          call FindOCIndex(iDOC, &
                           iSOC, &
                           SGS%CDCInfo(nOpen,mult)%OCDOC(:,iRefSm), &
                           SGS%CDCInfo(nOpen,mult)%OCSOC(:,iRefSm), &
                           nOC, &
                           iOC)

          if (SG_Checks) then  ! sanity check, only when verifying CI developments
            if (SGS%CDCInfo(nOpen,Mult)%OCReverse(iDOC,iSOC,iRefSm) /= iOC) then
              write(u6,'(A)') 'SGS%CDCInfo(nOpen,Mult)%OCReverse(iDOC,iSOC,iRefSm) /= iOC'
              call Abend()
            end if
          end if

          call FindSpinIndex(OpenSP, &
                             iOpen, &
                             SGS%CDCInfo(iOpen,mult)%DetSpin, &
                             iSpinPattern)

          if (SG_Checks) then  ! sanity check, only when verifying CI developments
            do i=1,iOpen
              if (OpenSP(i) == 1_iwp) then
                if (SGS%CDCInfo(iOpen,mult)%DetSpin(i,iSpinPattern) /= 1_iwp) then
                  write(u6,'(A)') 'Mk_SDREO_SGUGA_New: alpha-spin mismatch'
                  call Abend()
                end if
              else
                if (SGS%CDCInfo(iOpen,mult)%DetSpin(i,iSpinPattern) /= 0_iwp) then
                  write(u6,'(A)') 'Mk_SDREO_SGUGA_New: beta-spin mismatch'
                  call Abend()
                end if
              end if
            end do
          end if
          iAdrSDConfOrder = iOpenOffset(iOpen) &
                            +(iOC-1_iwp)*SGS%CDCInfo(iOpen,mult)%nDet+iSpinPattern
          iDet = iDet+1_iwp
          if (SG_Checks) then  ! sanity check, only when verifying CI developments
            if (iAdrSDConfOrder < 1_iwp .or. &
                iAdrSDConfOrder > nDetTotal) then
              write(u6,'(A,2(1X,I0))') 'ADDRESS RANGE',iAdrSDConfOrder,nDetTotal
              write(u6,'(A)') 'Mk_SDREO_SGUGA_New: invalid configuration address'
              call Abend()
            end if
            if (SGS%SDREO(iRefSm)%List(iAdrSDConfOrder) /= 0_iwp) then
              write(u6,'(A,3(1X,I0))') 'DUPLICATE SDREO',iRefSm,iAdrSDConfOrder,iDet
              write(u6,'(A)') 'Mk_SDREO_SGUGA: duplicate SDREO address'
              call Abend()
            end if
          end if

          SGS%SDREO(iRefSm)%List(iAdrSDConfOrder) = iSign*iDet

        end do
      end do

      call mma_deallocate(AlphaStrings)
      call mma_deallocate(BetaStrings)

    end do
    if (SG_Checks) then  ! sanity check, only when verifying CI developments
      if (iDet /= nDetTotal) then
        write(u6,'(A,3(1X,I0))') 'SDREO SIZE MISMATCH',iRefSm,iDet,nDetTotal
        write(u6,'(A)') 'Mk_SDREO_SGUGA: incorrect determinant count'
        call Abend()
      end if

      if (any(SGS%SDREO(iRefSm)%List == 0_iwp)) then
        write(u6,'(A,1X,I0)') 'INCOMPLETE SDREO FOR SYMMETRY',iRefSm
        write(u6,'(A)') 'Mk_SDREO_SGUGA: incomplete SDREO vector'
        call Abend()
      end if
    end if

  end do

contains

  subroutine FindSpinIndex(OpenSP,iOpen,DetSpin,iSpinPattern)

    implicit none

    integer(kind=iwp), intent(In) :: iOpen
    integer(kind=iwp), intent(In) :: OpenSP(max(1_iwp,iOpen))
    integer(kind=iwp), intent(In) :: DetSpin(:,:)
    integer(kind=iwp), intent(Out) :: iSpinPattern

    integer(kind=iwp) :: i
    integer(kind=iwp) :: j
    integer(kind=iwp) :: nSpin
    integer(kind=iwp) :: nMatch
    integer(kind=iwp) :: SpinValue
    logical :: Match

    nSpin = size(DetSpin,2)
    iSpinPattern = 0_iwp
    nMatch = 0_iwp

    if (iOpen == 0_iwp) then

      if (nSpin /= 1_iwp) then
        write(u6,'(A)') 'FindSpinIndex: invalid closed-shell spin table'
        call Abend()
      end if

      iSpinPattern = 1_iwp
      return

    end if

    if (size(DetSpin,1) < iOpen) then
      write(u6,'(A)') 'FindSpinIndex: spin table has too few rows'
      call Abend()
    end if

    SpinValue = -1_iwp
    do j=1,nSpin

      Match = .true.

      do i=1,iOpen

        if (OpenSP(i) == 1_iwp) then
          SpinValue = 1_iwp
        else if (OpenSP(i) == -1_iwp) then
          SpinValue = 0_iwp
        else
          write(u6,'(A)') 'FindSpinIndex: invalid open-shell spin'
          call Abend()
        end if

        if (DetSpin(i,j) /= SpinValue) then
          Match = .false.
          exit
        end if

      end do

      if (Match) then
        iSpinPattern = j
        nMatch = nMatch+1_iwp
      end if

    end do

    if (nMatch == 0_iwp) then
      write(u6,'(A)') 'FindSpinIndex: open-spin pattern not found'
      write(u6,'(A,2(1X,I0))') 'IOPEN NSPIN',iOpen,nSpin
      write(u6,'(A)') 'FindSpinIndex: spin pattern not found'
      call Abend()
    end if

    if (nMatch > 1_iwp) then
      write(u6,'(A,1X,I0)') 'FindSpinIndex: duplicate count',nMatch
      write(u6,'(A)') 'FindSpinIndex: duplicate spin pattern'
      call Abend()
    end if

  end subroutine FindSpinIndex

  subroutine FindOCIndex(iDOC,iSOC,OCDOC,OCSOC,nOC,iOC)

    implicit none

    integer(kind=iwp), intent(In) :: iDOC
    integer(kind=iwp), intent(In) :: iSOC
    integer(kind=iwp), intent(In) :: nOC
    integer(kind=iwp), intent(In) :: OCDOC(:)
    integer(kind=iwp), intent(In) :: OCSOC(:)
    integer(kind=iwp), intent(Out) :: iOC

    integer(kind=iwp) :: i
    integer(kind=iwp) :: nMatch

    if (nOC < 1_iwp) then
      write(u6,'(A)') 'FindOCIndex: no configurations for symmetry'
      call Abend()
    end if

    if (nOC > size(OCDOC) .or. nOC > size(OCSOC)) then
      write(u6,'(A)') 'FindOCIndex: invalid configuration count'
      call Abend()
    end if

    iOC = 0_iwp
    nMatch = 0_iwp

    do i=1,nOC

      if (OCDOC(i) /= iDOC) cycle
      if (OCSOC(i) /= iSOC) cycle

      iOC = i
      nMatch = nMatch+1_iwp

    end do

    if (nMatch == 0_iwp) then
      write(u6,'(A)') 'FindOCIndex: occupation configuration not found'
      write(u6,'(A,2(1X,I0))') 'IDOC ISOC',iDOC,iSOC
      write(u6,'(A,1X,I0)') 'NOC',nOC
      write(u6,'(A)') 'FindOCIndex: configuration not found'
      call Abend()
    end if

    if (nMatch > 1_iwp) then
      write(u6,'(A)') 'FindOCIndex: duplicate occupation configuration'
      write(u6,'(A,2(1X,I0))') 'IDOC ISOC',iDOC,iSOC
      write(u6,'(A,1X,I0)') 'NMATCH',nMatch
      write(u6,'(A)') 'FindOCIndex: duplicate configuration'
      call Abend()
    end if

  end subroutine FindOCIndex

  subroutine FindSOCIndex(iDetOC,nElec,SOCTable,iSOC)

    implicit none

    integer(kind=iwp), intent(In) :: nElec
    integer(kind=iwp), intent(In) :: iDetOC(nElec)
    integer(kind=iwp), intent(In) :: SOCTable(:,:)
    integer(kind=iwp), intent(Out) :: iSOC

    integer(kind=iwp) :: i
    integer(kind=iwp) :: j
    integer(kind=iwp) :: nSOC
    logical :: Match

    nSOC = size(SOCTable,1)
    iSOC = 0_iwp

    do j=1,size(SOCTable,2)

      Match = .true.

      do i=1,nSOC
        if (count(iDetOC == SOCTable(i,j)) /= 1) then
          Match = .false.
          exit
        end if
      end do

      if (Match) then
        iSOC = j
        exit
      end if

    end do

    if (iSOC == 0_iwp) then
      Error stop 'FindSOCIndex: SOC list not found'
    end if

  end subroutine FindSOCIndex

  subroutine FindDOCIndex(iDetOC,nElec,DOCTable,iDOC)

    use Definitions, only: iwp

    implicit none

    integer(kind=iwp), intent(In) :: nElec
    integer(kind=iwp), intent(In) :: iDetOC(nElec)
    integer(kind=iwp), intent(In) :: DOCTable(:,:)
    integer(kind=iwp), intent(Out) :: iDOC

    integer(kind=iwp) :: DOCOcc(size(DOCTable,1))
    integer(kind=iwp) :: nDOC
    integer(kind=iwp) :: nDOCTable
    integer(kind=iwp) :: iElec
    integer(kind=iwp) :: jDOC

    nDOC = 0_iwp
    nDOCTable = size(DOCTable,2)
    iElec = 1_iwp
    iDOC = 0_iwp

    do while (iElec <= nElec)

      if (iElec < nElec) then
        if (iDetOC(iElec) == iDetOC(iElec+1_iwp)) then

          nDOC = nDOC+1_iwp

          if (nDOC > size(DOCOcc)) then
            Error stop 'FindDOCIndex: too many double occupations'
          end if

          DOCOcc(nDOC) = iDetOC(iElec)
          iElec = iElec+2_iwp
          cycle

        end if
      end if

      iElec = iElec+1_iwp

    end do

    if (nDOC == 0_iwp) then

      if (size(DOCTable,1) /= 0) then
        Error stop 'FindDOCIndex: expected zero-length DOC table'
      end if

      if (size(DOCTable,2) /= 1) then
        Error stop 'FindDOCIndex: invalid empty DOC table'
      end if

      iDOC = 1_iwp
      return

    end if

    if (nDOC /= size(DOCTable,1)) then
      Error stop 'FindDOCIndex: inconsistent number of double occupations'
    end if

    do jDOC=1,nDOCTable
      if (all(DOCOcc(1:nDOC) == DOCTable(:,jDOC))) then
        iDOC = jDOC
        exit
      end if
    end do

    if (iDOC == 0_iwp) then
      Error stop 'FindDOCIndex: DOC occupation not found'
    end if

  end subroutine FindDOCIndex

  subroutine ExtractOccupationPattern(iDetOC,iDetSP,nElec, &
                                      DOCOcc,SOCOcc,OpenSP, &
                                      nDouble,nOpen)

    use Definitions, only: iwp

    implicit none

    integer(kind=iwp), intent(In) :: nElec
    integer(kind=iwp), intent(In) :: iDetOC(nElec)
    integer(kind=iwp), intent(In) :: iDetSP(nElec)

    integer(kind=iwp), intent(Out) :: DOCOcc(nElec)
    integer(kind=iwp), intent(Out) :: SOCOcc(nElec)
    integer(kind=iwp), intent(Out) :: OpenSP(nElec)
    integer(kind=iwp), intent(Out) :: nDouble
    integer(kind=iwp), intent(Out) :: nOpen

    integer(kind=iwp) :: iElec

    DOCOcc = 0_iwp
    SOCOcc = 0_iwp
    OpenSP = 0_iwp

    nDouble = 0_iwp
    nOpen = 0_iwp

    iElec = 1_iwp

    do while (iElec <= nElec)

      if (iElec < nElec) then

        if (iDetOC(iElec) == iDetOC(iElec+1_iwp)) then

          ! A repeated orbital number represents a doubly
          ! occupied orbital. BuildOrderedDeterminant places
          ! the alpha electron before the beta electron.

          if (iDetSP(iElec) /= 1_iwp .or. &
              iDetSP(iElec+1_iwp) /= -1_iwp) then
            Error stop &
              'ExtractOccupationPattern: invalid double occupation'
          end if

          nDouble = nDouble+1_iwp
          DOCOcc(nDouble) = iDetOC(iElec)

          iElec = iElec+2_iwp
          cycle

        end if

      end if

      ! A non-repeated orbital number represents a singly
      ! occupied orbital.

      if (abs(iDetSP(iElec)) /= 1_iwp) then
        Error stop 'ExtractOccupationPattern: invalid spin label'
      end if

      nOpen = nOpen+1_iwp
      SOCOcc(nOpen) = iDetOC(iElec)
      OpenSP(nOpen) = iDetSP(iElec)

      iElec = iElec+1_iwp

    end do

! Every electron must belong either to a doubly occupied
! orbital or to a singly occupied orbital.

    if (2_iwp*nDouble+nOpen /= nElec) then
      Error stop 'ExtractOccupationPattern:*inconsistent electron count'
    end if

  end subroutine ExtractOccupationpattern

  subroutine BuildOrderedDeterminant(AlphaString,BetaString, &
                                     nAlpha,nBeta, &
                                     iDetOC,iDetSP,iSign)

    use Definitions, only: iwp

    implicit none

    integer(kind=iwp), intent(In) :: nAlpha
    integer(kind=iwp), intent(In) :: nBeta
    integer(kind=iwp), intent(In) :: AlphaString(max(1_iwp,nAlpha))
    integer(kind=iwp), intent(In) :: BetaString(max(1_iwp,nBeta))

    integer(kind=iwp), intent(Out) :: iDetOC(nAlpha+nBeta)
    integer(kind=iwp), intent(Out) :: iDetSP(nAlpha+nBeta)
    integer(kind=iwp), intent(Out) :: iSign

    integer(kind=iwp) :: iAlpha
    integer(kind=iwp) :: iBeta
    integer(kind=iwp) :: iElec
    integer(kind=iwp) :: nCross

    iAlpha = 1_iwp
    iBeta = 1_iwp
    iSign = 1_iwp

    do iElec=1,nAlpha+nBeta

      if (iAlpha <= nAlpha .and. iBeta <= nBeta) then

        if (AlphaString(iAlpha) <= BetaString(iBeta)) then

          ! Add the next alpha electron.  The use of <= ensures
          ! that alpha precedes beta in a doubly occupied orbital.

          iDetOC(iElec) = AlphaString(iAlpha)
          iDetSP(iElec) = 1_iwp

          iAlpha = iAlpha+1_iwp

        else

          ! Add the next beta electron.

          iDetOC(iElec) = BetaString(iBeta)
          iDetSP(iElec) = -1_iwp

          ! This beta electron passes all remaining alpha electrons.

          nCross = nAlpha-iAlpha+1_iwp

          if (mod(nCross,2_iwp) /= 0_iwp) then
            iSign = -iSign
          end if

          iBeta = iBeta+1_iwp

        end if

      else if (iBeta > nBeta) then

        ! No beta electrons remain.

        iDetOC(iElec) = AlphaString(iAlpha)
        iDetSP(iElec) = 1_iwp

        iAlpha = iAlpha+1_iwp

      else

        ! No alpha electrons remain.

        iDetOC(iElec) = BetaString(iBeta)
        iDetSP(iElec) = -1_iwp

        iBeta = iBeta+1_iwp

      end if

    end do

! Basic consistency checks.

    do iElec=2,nAlpha+nBeta
      if (iDetOC(iElec) < iDetOC(iElec-1_iwp)) then
        Error stop 'BuildOrderedDeterminant: unordered occupation list'
      end if

      if (iDetOC(iElec) == iDetOC(iElec-1_iwp)) then
        if (iDetSP(iElec-1_iwp) /= 1_iwp .or. &
            iDetSP(iElec) /= -1_iwp) then
          Error stop 'BuildOrderedDeterminant: invalid double occupation'
        end if
      end if
    end do

  end subroutine BuildOrderedDeterminant

  subroutine Count_Strings_By_Symmetry(SGS,nElec,nStringSym)

    implicit none

    type(SGStruct), intent(in) :: SGS

    integer(kind=iwp), intent(in) :: nElec
    integer(kind=iwp), intent(out) :: nStringSym(SGS%nSym)

    integer(kind=iwp) :: iOrb
    integer(kind=iwp) :: iSym
    integer(kind=iwp) :: nOcc
    integer(kind=iwp) :: newSym

    integer(kind=iwp) :: nString(0:nElec,SGS%nSym)
    integer(kind=iwp) :: nStringNew(0:nElec,SGS%nSym)

    ! nString(nOcc,iSym) is the number of strings constructed
    ! from the orbitals processed so far that contain nOcc
    ! electrons and have spatial symmetry iSym.
    nString(:,:) = 0_iwp

    ! The empty string is totally symmetric.
    nString(0,1) = 1_iwp

    do iOrb=1,SGS%nLev

      ! Contributions in which orbital iOrb is unoccupied.
      nStringNew(:,:) = nString(:,:)

      ! Contributions in which orbital iOrb is occupied.
      do nOcc=0,min(nElec-1_iwp,iOrb-1_iwp)
        do iSym=1,SGS%nSym

          if (nString(nOcc,iSym) == 0_iwp) cycle

          newSym = Mul(iSym,SGS%ISm(iOrb))

          nStringNew(nOcc+1,newSym) = &
            nStringNew(nOcc+1,newSym)+ &
            nString(nOcc,iSym)

        end do
      end do

      nString(:,:) = nStringNew(:,:)

    end do

    nStringSym(:) = nString(nElec,:)

  end subroutine Count_Strings_By_Symmetry

  subroutine GenerateStringSet(nElec,targetSym,nStrings,Strings)

    integer(kind=iwp), intent(in) :: nElec
    integer(kind=iwp), intent(in) :: targetSym
    integer(kind=iwp), intent(in) :: nStrings

    integer(kind=iwp), intent(out) :: &
      Strings(max(1_iwp,nElec),nStrings)

    integer(kind=iwp), allocatable :: String(:)

    integer(kind=iwp) :: iString
    integer(kind=iwp) :: nOrb

    nOrb = SGS%nLev

    Strings = 0_iwp
    iString = 0_iwp

    ! The empty string has totally symmetric spatial symmetry.
    if (nElec == 0_iwp) then

      if (targetSym == 1_iwp) then
        if (nStrings /= 1_iwp) then
          error stop 'GenerateStringSet: invalid empty-string count'
        end if
      else
        if (nStrings /= 0_iwp) then
          error stop 'GenerateStringSet: invalid empty-string symmetry'
        end if
      end if

      return

    end if

    if (nElec < 0_iwp .or. nElec > nOrb) then
      error stop 'GenerateStringSet: invalid electron count'
    end if

    call mma_allocate(String,nElec,Label='String')
    String = 0_iwp

    ! Generate combinations in the same colexicographic order
    ! used by LUCIA.
    call GenerateLevel(nElec,nOrb,nElec,targetSym,nStrings, &
                       String,Strings,iString)

    if (iString /= nStrings) then
      error stop 'GenerateStringSet: incorrect number of strings'
    end if

    call mma_deallocate(String)

  end subroutine GenerateStringSet

  recursive subroutine GenerateLevel(nChoose,maxOrb,nElec,targetSym,nStrings, &
                                     String,Strings,iString)

    integer(kind=iwp), intent(in) :: nChoose
    integer(kind=iwp), intent(in) :: maxOrb
    integer(kind=iwp), intent(in) :: nElec
    integer(kind=iwp), intent(in) :: targetSym
    integer(kind=iwp), intent(in) :: nStrings

    integer(kind=iwp), intent(inout) :: String(nElec)
    integer(kind=iwp), intent(out) :: Strings(nElec,nStrings)
    integer(kind=iwp), intent(inout) :: iString

    integer(kind=iwp) :: i
    integer(kind=iwp) :: iOrb
    integer(kind=iwp) :: iSym

    if (nChoose == 0_iwp) then

      iSym = 1_iwp

      do i=1,nElec
        iSym = Mul(iSym,SGS%ISm(String(i)))
      end do

      if (iSym == targetSym) then

        iString = iString+1_iwp

        if (iString > nStrings) then
          error stop 'GenerateLevel: too many strings'
        end if

        Strings(:,iString) = String(:)

      end if

      return

    end if

    ! Generate strings in colexicographic order. The highest
    ! occupied orbital is the outermost ordering criterion.
    do iOrb=nChoose,maxOrb

      String(nChoose) = iOrb

      call GenerateLevel( &
        nChoose-1_iwp,iOrb-1_iwp, &
        nElec,targetSym,nStrings, &
        String,Strings,iString)

    end do

  end subroutine GenerateLevel

end subroutine Mk_SDREO_SGUGA

pure function SOC_Sym(ISm,nSOC,SOCOcc) result(iSym)
  use Symmetry_Info, only: Mul
  use Definitions, only: iwp
  implicit none

  integer(iwp), intent(in) :: ISm(:)
  integer(iwp), intent(in) :: nSOC
  integer(iwp), intent(in) :: SOCOcc(:)

  integer(iwp) :: iSym
  integer(iwp) :: i

  iSym = 1
  do i=1,nSOC
    iSym = Mul(iSym,ISm(SOCOcc(i)))
  end do

end function SOC_Sym

subroutine Count_OpenShell_OCPs(ISm,nLev,nOpen,nSym,nOpenOCP)

  use Symmetry_Info, only: Mul
  use Definitions, only: iwp

  implicit none

  integer(iwp), intent(in) :: nLev
  integer(iwp), intent(in) :: nOpen
  integer(iwp), intent(in) :: nSym
  integer(iwp), intent(in) :: ISm(nLev)

  integer(iwp), intent(out) :: nOpenOCP(nSym)

  nOpenOCP(:) = 0

  if (nOpen == 0) then
    nOpenOCP(1) = 1
    return
  end if

  call recurse(1,1,1)

  if (sum(nOpenOCP) /= binomial(nOpen,nLev)) then
    write(u6,*) 'Count_OpenShell_OCPs failed'
    write(u6,*) 'NOPEN = ',nOpen
    write(u6,*) 'SUM(NOPENOCP) = ',sum(nOpenOCP)
    write(u6,*) 'EXPECTED      = ',binomial(nOpen,nLev)
    call Abend()
  end if

  if (any(nOpenOCP < 0)) then
    write(u6,*) 'Count_OpenShell_OCPs failed'
    write(u6,*) 'Negative count encountered'
    call Abend()
  end if

contains

  recursive subroutine recurse(iPicked,iStart,iSymProd)

    integer(iwp), intent(in) :: iPicked
    integer(iwp), intent(in) :: iStart
    integer(iwp), intent(in) :: iSymProd

    integer(iwp) :: i

    if (iPicked > nOpen) then
      nOpenOCP(iSymProd) = nOpenOCP(iSymProd)+1
      return
    end if

    do i=iStart,nLev-(nOpen-iPicked)

      call recurse( &
        iPicked+1, &
        i+1, &
        Mul(iSymProd,ISm(i)))

    end do

  end subroutine recurse

end subroutine Count_OpenShell_OCPs

subroutine sg_csfdet(nopen, &
                     ndet, &
                     detspin, &
                     ncsf, &
                     twos, &
                     cdc)

  use Constants, only: Half, One, Two, Zero
  use Definitions, only: wp, iwp

  implicit none

  integer(iwp), intent(in) :: nopen
  integer(iwp), intent(in) :: ndet
  integer(iwp), intent(in) :: ncsf
  integer(iwp), intent(in) :: detspin(nopen,ndet)
  integer(iwp), intent(in) :: twos(0:nopen,ncsf)

  real(wp), intent(out) :: cdc(ndet,ncsf)

  integer(iwp) :: idet
  integer(iwp) :: icsf
  integer(iwp) :: i

  real(wp) :: coef, norm
  real(wp) :: sign
  real(wp) :: s
  real(wp) :: m(max(1,nopen))

  logical :: up

#ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) '========================================'
  write(u6,*) ' SG_CSFDET'
  write(u6,*) '========================================'
  write(u6,*) ' NOPEN, nDet, nCSF = ',nopen,ndet,ncsf
  write(u6,*)
#endif

  do idet=1,ndet
    call det_mchain(nopen,detspin(:,idet),m)

#ifdef _DEBUGPRINT_
    write(u6,*)
    write(u6,*) '----------------------------------------'
    write(u6,*) ' Determinant ',idet
    write(u6,*) '----------------------------------------'
    write(u6,'(A,100I3)') ' DETSPIN = ',(detspin(i,idet),i=1,nopen)
    write(u6,'(A,100F8.3)') ' MCHAIN  = ',(m(i),i=1,nopen)
#endif

    do icsf=1,ncsf

#ifdef _DEBUGPRINT_
      write(u6,*)
      write(u6,*) ' CSF = ',icsf
      write(u6,'(A,100I3)') ' TWOS = ',(twos(i,icsf),i=0,nopen)
      write(u6,'(A)',advance='no')
      write(u6,*) ' SCHAIN:'
      do i=1,nopen
        write(u6,'(I4,F10.4)') i,0.5_wp*real(twos(i,icsf),wp)
      end do
#endif

      coef = One
      sign = One

      do i=1,nopen
        s = Half*real(twos(i,icsf),wp)
        ! walk step:
        ! twos(i)=twos(i-1)+1  => up
        ! twos(i)=twos(i-1)-1  => down

        up = (twos(i,icsf) > twos(i-1,icsf))
        if (up) then
          if (detspin(i,idet) == 1) then
            coef = coef*(s+m(i))/(Two*s)
          else
            coef = coef*(s-m(i))/(Two*s)
          end if

        else

          if (detspin(i,idet) == 1) then
            coef = coef*(s-m(i)+One)/(Two*s+Two)
            sign = -sign
          else
            coef = coef*(s+m(i)+One)/(Two*s+Two)
          end if

        end if

#ifdef _DEBUGPRINT_
        write(u6,*)
        write(u6,*) ' Step i = ',i
        write(u6,*) ' twos(i-1) = ',twos(i-1,icsf)
        write(u6,*) ' twos(i)   = ',twos(i,icsf)
        write(u6,*) ' S = ',s
        write(u6,*) ' M = ',m(i)

        if (up) then
          write(u6,*) ' walk step = UP'
        else
          write(u6,*) ' walk step = DOWN'
        end if

        if (detspin(i,idet) == 1) then
          write(u6,*) ' determinant spin = ALPHA'
        else
          write(u6,*) ' determinant spin = BETA'
        end if
#endif

      end do

#ifdef _DEBUGPRINT_
      write(u6,*)
      write(u6,*) ' FINAL PRODUCT = ',coef
      write(u6,*) ' FINAL SIGN    = ',sign
#endif

      if (coef < Zero) then
        cdc(idet,icsf) = Zero
      else
        cdc(idet,icsf) = sign*sqrt(coef)
      end if

#ifdef _DEBUGPRINT_
      write(u6,*) ' CDC(',idet,',',icsf,') = ',cdc(idet,icsf)
#endif

    end do
  end do

#ifdef _DEBUGPRINT_
  call RecPrt('Final CDC matrix',' ',cdc,ndet,ncsf)
#endif

  do icsf=1,ncsf
    norm = Zero

    do idet=1,ndet
      norm = norm+cdc(idet,icsf)**2
    end do

    if (abs(One-Norm) > 1.0e-12_wp) then
      write(u6,*)
      write(u6,*) ' CSF ',icsf,' norm = ',norm
      call Abend()
    end if

  end do

contains
  subroutine det_mchain(nopen,det,m)

    use Constants, only: Half
    use Definitions, only: wp, iwp

    implicit none

    integer(iwp), intent(in) :: nopen
    integer(iwp), intent(in) :: det(max(1,nopen))
    real(wp), intent(out) :: m(max(1,nopen))

    integer(iwp) :: i

    m(1) = real(det(1),wp)-Half

    do i=2,nopen
      m(i) = m(i-1)+real(det(i),wp)-Half
    end do

  end subroutine det_mchain

end subroutine sg_csfdet

#ifdef _TEST_CSFDET_
subroutine test_sg_csfdet()

  use Definitions, only: wp, iwp

  implicit none

  integer(iwp), parameter :: maxopen = 8
  integer(iwp) :: nopen
  integer(iwp) :: mult
  integer(iwp) :: ms2

  integer(iwp), allocatable :: detspin(:,:)
  integer(iwp), allocatable :: twos(:,:)

  real(wp), allocatable :: cdc(:,:)
  integer(iwp) :: ndet, ncsf

#ifdef _DEBUGPRINT_
  integer(iwp) :: idet, icsf

  character(len=1) :: walk(Maxopen)
  character(len=1) :: str(Maxopen)
  integer(iwp) :: i

  write(u6,*)
  write(u6,*) '===================================='
  write(u6,*) ' TEST SG_CSFDET'
  write(u6,*) '===================================='
  write(u6,*)
#endif

  do nopen=1,maxopen

#ifdef _DEBUGPRINT_
    write(u6,*)
    write(u6,*) '========================================='
    write(u6,*) 'NOPEN = ',nopen
    write(u6,*) '========================================='
#endif

    do mult=mod(nopen,2)+1,nopen+1,2

      ms2 = mult-1

#ifdef _DEBUGPRINT_
      write(u6,*)
      write(u6,*) 'MULT = ',mult
      write(u6,*) 'MS2  = ',ms2
#endif

      call generate_detspin(nopen,ms2,detspin,ndet)
#ifdef _DEBUGPRINT_
      write(u6,*)
      write(u6,*) 'Determinant strings'
      do idet=1,ndet
        do i=1,nopen
          if (detspin(i,idet) == 1) then
            str(i) = 'a'
          else
            str(i) = 'b'
          end if
        end do
        write(u6,'(I4,2X,100A1)') idet,(str(i),i=1,nopen)
      end do
#endif

      call generate_twos(nopen,mult,twos,ncsf)
#ifdef _DEBUGPRINT_
      write(u6,*)
      write(u6,*) 'CSF walks'
      do icsf=1,ncsf
        do i=1,nopen
          if (twos(i,icsf) > twos(i-1,icsf)) then
            walk(i) = 'u'
          else
            walk(i) = 'd'
          end if
        end do
        write(u6,'(I4,2X,100A1)') icsf,(walk(i),i=1,nopen)
      end do
#endif

#ifdef _TEST_CODE_
      if (ms2 == mult-1) then
        call check_twos(nopen,ndet,detspin,ncsf,twos)
      end if
#endif

      call mma_allocate(cdc,nDet,nCSF,Label='cdc')
      call sg_csfdet(nopen,ndet,detspin,ncsf,twos,cdc)

      call mma_deallocate(detspin)
      call mma_deallocate(twos)
      call mma_deallocate(cdc)

    end do

  end do

end subroutine test_sg_csfdet
#endif

subroutine generate_detspin(nopen,ms2,detspin,ndet)

  use Definitions, only: iwp
  use stdalloc, only: mma_allocate

  implicit none

  integer(iwp), intent(in) :: nopen
  integer(iwp), intent(in) :: ms2

  integer(iwp), intent(out) :: ndet
  integer(iwp), allocatable, intent(out) :: detspin(:,:)

  integer(iwp) :: ibit
  integer(iwp) :: idet
  integer(iwp) :: i
  integer(iwp) :: mx

!
! Count determinants
!

  ndet = 0
  mx = 2**nopen

  do ibit=0,mx-1

    if (2*popcnt(ibit)-nopen == ms2) then
      ndet = ndet+1
    end if

  end do

  call mma_allocate(detspin,[1,max(1,nopen)],[1,ndet],Label='detspin')

!
! Build determinant list
!

  idet = 0

  do ibit=0,mx-1

    if (2*popcnt(ibit)-nopen /= ms2) cycle

    idet = idet+1

    do i=1,nopen
      detspin(i,idet) = ibits(ibit,i-1,1)
    end do

  end do

end subroutine generate_detspin

subroutine generate_twos(nopen,mult,twos,ncsf)

  use Definitions, only: iwp
  use stdalloc, only: mma_allocate

  implicit none

  integer(iwp), intent(in) :: nopen
  integer(iwp), intent(in) :: mult

  integer(iwp), intent(out) :: ncsf
  integer(iwp), allocatable, intent(out) :: twos(:,:)

  integer(iwp) :: ibit
  integer(iwp) :: i
  integer(iwp) :: ms2l
  integer(iwp) :: mx
  integer(iwp) :: icsf

  logical :: upper

!
! Count walks
!

  ncsf = 0
  mx = 2**nopen

  do ibit=0,mx-1
    ms2l = 0
    upper = .true.
    do i=1,nopen
      if (ibits(ibit,nopen-i,1) == 1) then
        ms2l = ms2l+1
      else
        ms2l = ms2l-1
      end if

      if (ms2l < 0) then
        upper = .false.
        exit
      end if
    end do

    if (.not. upper) cycle
    if (ms2l /= mult-1) cycle
    ncsf = ncsf+1

  end do

  call mma_allocate(twos,[0,nopen],[1,ncsf],Label='twos')

!
! Build walks
!

  icsf = 0
  do ibit=0,mx-1

    ms2l = 0
    upper = .true.
    do i=1,nopen
      if (ibits(ibit,i-1,1) == 1) then
        ms2l = ms2l+1
      else
        ms2l = ms2l-1
      end if

      if (ms2l < 0) then
        upper = .false.
        exit
      end if
    end do

    if (.not. upper) cycle
    if (ms2l /= mult-1) cycle

    icsf = icsf+1
    twos(0,icsf) = 0
    ms2l = 0

    do i=1,nopen
      if (ibits(ibit,i-1,1) == 1) then
        ms2l = ms2l+1
      else
        ms2l = ms2l-1
      end if

      twos(i,icsf) = ms2l
    end do
  end do

end subroutine generate_twos

#ifdef _TEST_CODE_
subroutine check_twos(nopen,ndet,detspin,ncsf,twos)

  use Definitions, only: iwp, u6

  implicit none

  integer(iwp), intent(in) :: nopen
  integer(iwp), intent(in) :: ndet
  integer(iwp), intent(in) :: ncsf

  integer(iwp), intent(in) :: detspin(nopen,ndet)
  integer(iwp), intent(in) :: twos(0:nopen,ncsf)

  integer(iwp) :: idet
  integer(iwp) :: iup
  integer(iwp) :: i
  integer(iwp) :: ms2l
  integer(iwp) :: step

  logical :: upper

  iup = 0

  do idet=1,ndet

    upper = .true.
    ms2l = 0

    do i=1,nopen

      if (detspin(i,idet) == 1) then
        ms2l = ms2l+1
      else
        ms2l = ms2l-1
      end if

      if (ms2l < 0) then
        upper = .false.
        exit
      end if

    end do

    if (.not. upper) cycle

    iup = iup+1

    if (iup > ncsf) then
      write(u6,*) 'CHECK_TWOS FAILED'
      write(u6,*) 'Too many upper determinants'
      write(u6,*) 'IUP  = ',iup
      write(u6,*) 'NCSF = ',ncsf
      call Abend()
    end if

    do i=1,nopen

      step = twos(i,iup)-twos(i-1,iup)

      if (detspin(i,idet) == 1) then

        if (step /= 1) then
          write(u6,*) 'CHECK_TWOS FAILED'
          write(u6,*) 'IDET = ',idet
          write(u6,*) 'IUP  = ',iup
          write(u6,*) 'POSITION = ',i
          call Abend()
        end if

      else

        if (step /= -1) then
          write(u6,*) 'CHECK_TWOS FAILED'
          write(u6,*) 'IDET = ',idet
          write(u6,*) 'IUP  = ',iup
          write(u6,*) 'POSITION = ',i
          call Abend()
        end if

      end if

    end do

  end do

  if (iup /= ncsf) then
    write(u6,*) 'CHECK_TWOS FAILED'
    write(u6,*) 'Number of upper determinants /= number of CSFs'
    write(u6,*) 'IUP  = ',iup
    write(u6,*) 'NCSF = ',ncsf
    call Abend()
  end if

end subroutine check_twos
#endif

#ifdef _GOSPER_
!Todo: implement in generate_detspin and generate_twos
subroutine next_combination(x)
  implicit none
  integer, intent(inout) :: x

  integer :: c, r

  c = iand(x,-x)
  r = x+c
  x = ior(r,ishft(ieor(r,x)/c,-2))

end subroutine next_combination
integer :: x, m, n

m = 5
n = 3

x = 2**n-1

do while (x < 2**m)

  write(u6,'(I4,2X,B5.5)') x,x

  call next_combination(x)

end do
#endif

!***********************************************************************
! Generate the next orbital occupation configuration in Lucia order.
!
! Occupations are represented by the occupied orbital index of each
! electron:
!
!    OC = (1,1,2,2,3)
!
! corresponds to
!
!    DOC = (1,2)
!    SOC = (3)
!
! The occupation list is stored in nondecreasing order and no orbital
! may occur more than twice. With INI=1 the first configuration is
! generated. With INI=0 the next configuration in lexical Lucia order
! is generated.
!
! Input:
!    nElec  - number of electrons
!    nLev   - number of active orbitals
!    Ini    - 1: initialize, 0: increment
!
! Output:
!    OC     - orbital occupation configuration
!    NoNew  - 0: new configuration generated
!             1: no further configuration exists
!
!***********************************************************************
subroutine Nxt_OC(OC,nElec,nLev,Ini,NoNew)

  implicit none

  integer, intent(In) :: nElec
  integer, intent(In) :: nLev
  integer, intent(In) :: Ini
  integer, intent(InOut) :: OC(nElec)
  integer, intent(Out) :: NoNew

  integer :: I
  integer :: IAdd
  integer :: IEL
  integer :: Increase
  integer :: JOrb
  integer :: nDouble

  if (Ini == 1) then
    if (nElec <= 2*nLev) then
      NoNew = 0
      nDouble = nElec/2
      do I=1,nDouble
        OC(2*I-1) = I
        OC(2*I) = I
      end do
      if (2*nDouble /= nElec) &
        OC(2*nDouble+1) = nDouble+1
    else
      NoNew = 1
    end if
  else
    IAdd = 1
    IEL = 0
    do while (IAdd == 1)
      IEL = IEL+1
      Increase = 0
      if (IEL < nElec) then
        if (OC(IEL) < OC(IEL+1)-1) then
          Increase = 1
        else if (OC(IEL) == OC(IEL+1)-1) then
          if (IEL == nElec-1) then
            Increase = 1
          else if (OC(IEL+1) /= OC(IEL+2)) then
            Increase = 1
          end if
        end if
      else if (IEL == nElec) then
        if (OC(IEL) < nLev) then
          Increase = 1
        else
          NoNew = 1
          exit
        end if
      else
        NoNew = 1
        exit
      end if

      if (Increase == 1) then
        NoNew = 0
        OC(IEL) = OC(IEL)+1
        nDouble = (IEL-1)/2
        do JOrb=1,nDouble
          OC(2*JOrb-1) = JOrb
          OC(2*JOrb) = JOrb
        end do
        if (2*nDouble < IEL-1) &
          OC(IEL-1) = nDouble+1
        IAdd = 0
      end if
    end do
  end if

end subroutine Nxt_OC

subroutine Decode_OC(OC,nElec,nDOC,DOCOcc,nSOC,SOCOcc)

  implicit none

  integer, intent(In) :: nElec
  integer, intent(In) :: OC(nElec)

  integer, intent(Out) :: nDOC
  integer, intent(Out) :: nSOC
  integer, intent(Out) :: DOCOcc(:)
  integer, intent(Out) :: SOCOcc(:)

  integer :: i

  nDOC = 0
  nSOC = 0
  i = 1

  do while (i <= nElec)

    if (i < nElec) then
      if (OC(i) == OC(i+1)) then
        nDOC = nDOC+1
        DOCOcc(nDOC) = OC(i)
        i = i+2
      else
        nSOC = nSOC+1
        SOCOcc(nSOC) = OC(i)
        i = i+1
      end if
    else
      nSOC = nSOC+1
      SOCOcc(nSOC) = OC(i)
      i = i+1
    end if

  end do

end subroutine Decode_OC

end module SGUGA_Dets
