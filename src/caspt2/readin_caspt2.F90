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
! Copyright (C) Steven Vancoillie                                      *
!***********************************************************************

module InputData
  !SVC: this module contains a data structure to keep all input variables.
  ! Note that I use the standard 'allocate' here against the appropriate
  ! Molcas practices. The reason is that these are (i) very small, and
  ! (ii) there is need for allocating complex things such as derived
  ! types, which are not supported with stdalloc. Hence, the infraction.

  use definitions, only: wp,iwp

  implicit none

#include "compiler_features.h"

  type States
    Integer(kind=iwp), allocatable :: State(:)
  end type

  type InputTable
    ! TITL      one line with a descriptive name
    Character(len=128) :: Title = ' '
    ! FILE      file to read CAS/RAS reference from
    Character(len=128) :: File = 'JOBIPH'
    ! MULT      the number of states, followed by the ID of each state
    Logical(kind=iwp) :: MULT = .false.
    Integer(kind=iwp) :: nMultState = 0
    type(States)      :: MultGroup
    Logical(kind=iwp) :: AllMult = .false.
    ! XMUL      extended multi-state caspt2
    Logical(kind=iwp) :: XMUL = .false.
    Integer(kind=iwp) :: nXMulState = 0
    type(States)      :: XMulGroup
    Logical(kind=iwp) :: AllXMult = .false.
    ! RMUL      rotated multi-state caspt2
    Logical(kind=iwp) :: RMUL = .false.
    Integer(kind=iwp) :: nRMulState = 0
    type(States)      :: RMulGroup
    Logical(kind=iwp) :: AllRMult = .false.
    ! DWMS      use dynamical weighting to construct Fock
    Logical(kind=iwp) :: DWMS = .false.
    Integer(kind=iwp) :: DWType = 3
    Integer(kind=iwp) :: ZETA = 50
    ! EFOC      uses rotated E_0 energies with DWMS
    Logical(kind=iwp) :: EFOC = .false.
    ! LROO      compute only a single root, mutually exclusive with both MULT or XMUL
    Logical(kind=iwp) :: LROO = .false.
    Integer(kind=iwp) :: SingleRoot = 0
    ! RLXR      root for which the gradient is computed
    Integer(kind=iwp) :: RlxRoot = -1

    ! IPEA      sets the IP-EA shift
    Logical(kind=iwp) :: IPEA = .false.
    Real(kind=wp)     :: BSHIFT = 0.0_wp
    ! IMAG      size of extra 'imaginary' denominator shift
    Real(kind=wp)     :: ShiftI = 0.0_wp
    ! SHIF      size of extra denominator shift
    Real(kind=wp)     :: Shift = 0.0_wp

    ! several freeze-delete schemes, each of these should active
    ! the general flag below, to indicate additional conversion is
    ! needed on the input orbitals
    Logical(kind=iwp) :: modify_correlating_MOs = .false.
    ! AFRE      freeze orbitals that do not have sufficient density on specified 'active' atoms
    Logical(kind=iwp) :: aFreeze = .false.
    Integer(kind=iwp) :: lnFro = 0
    Real(kind=wp)     :: ThrFr = 0.0_wp,ThrDe = 0.0_wp
    Character(len=4),allocatable :: NamFro(:)
    ! LOVC      freeze orbitals that are not localized no the active site
    Logical(kind=iwp) :: LovCASPT2 = .false.
    Real(kind=wp)     :: Thr_Atm = 0.0_wp
    ! FNOC      delete a fraction of virtual orbitals
    Logical(kind=iwp) :: FnoCASPT2 = .false.
    Real(kind=wp)     :: VFrac = 0.0_wp
    ! DOMP
    Logical(kind=iwp) :: doMP2 = .false.
    ! DOEN
    Logical(kind=iwp) :: doEnv = .false.
    ! VIRA
    Logical(kind=iwp) :: VIRA = .false.
    ! GHOS      excludes ghost orbitals from the PT2 treatment
    Logical(kind=iwp) :: GhostDelete = .false.
    Real(kind=wp)     :: ThrGD = 0.0_wp

    ! FROZ      number of frozen orbitals in each irrep
    Logical(kind=iwp) :: FROZ = .false.
    Integer(kind=iwp),allocatable :: nFro(:)
    ! DELE      number of deleted orbitals in each irrep
    Logical(kind=iwp) :: DELE = .false.
    Integer(kind=iwp),allocatable :: nDel(:)
    ! DENS      computes full density matrix from the 1st-order wavefunction
    Logical(kind=iwp) :: DENS = .false.
    ! RFPE      make a perturbative reaction field calculation
    Logical(kind=iwp) :: RFPert = .false.
    ! THRE      thresholds for removal of:
    !   ThrsHN    zero-norm components in the first-order perturbed
    !             wave function
    !   ThrsHS    linear dependencies between components of the first-
    !             order perturbed wave function
    Logical(kind=iwp) :: THRE = .false.
    Real(kind=wp)     :: ThrsHN = 1.0e-10_wp,ThrsHS = 1.0e-8_wp
    ! MAXI      maximum number of iterations for solving a system of
    !           linear equations, default 20. A 0 indicates: use of
    !           the diagonal zeroth order hamiltonian
    Integer(kind=iwp) :: maxIter = 20
    ! Conv      convergence criteria for solving a system of linear equations
    Real(kind=wp)     :: ThrConv = 1.0e-6_wp
    ! NOMI      do not create an PM-CAS wavefunction file (JobMix)
    Logical(kind=iwp) :: NoMix = .false.
    ! NOMU      do not perform a multistate interaction
    Logical(kind=iwp) :: noMult = .false.
    ! ONLY      in a MS calculation, compute a single root with couplings to the other roots
    Integer(kind=iwp) :: OnlyRoot = 0
    ! EFFE      read Heff coupling terms from the input and perform only the multistate part
    Logical(kind=iwp) :: JMS = .false.
    Real(kind=wp),allocatable :: Heff(:,:)
    ! NOOR      do not print orbitals
    Logical(kind=iwp) :: PrOrb = .true.
    ! PROP      compute properties
    ! NOPR      do not compute properties
    Logical(kind=iwp) :: Properties = .false.
    ! transformation of reference (input) orbitals
    ! NOTR      do not transform to quasi-canonical orbitals,
    !           regardless of the state of the reference orbitals
    ! TRAN      transform to quasi-canonical orbitals, regardless
    !           of the state of the reference orbitals
    ! the default is to use transformation, unless the PT2 keyword
    ! was used in the rasscf program and the fock matrix is standard
    Character(len=8)  :: OrbIn = 'TRANSFOR'
    ! OFEM      add orbital-free embedding potential to hamiltonian
    Logical(kind=iwp) :: OFEmbedding = .false.
    ! OUTP      control extent of orbital printing
    Character(len=8)  :: OutFormat = 'DEFAULT '
    ! PRWF      print the CI coefficients above this threshold
    Real(kind=wp)     :: PrWF = 0.05_wp
    ! PRSD      print the determinant expansion of CSFs
    Logical(kind=iwp) :: PrSD = .false.
    ! NOOR      do not print any orbitals
    Logical(kind=iwp) :: NoOrb = .false.

    ! UNDOCUMENTED KEYWORDS
    ! CHOL
    Logical(kind=iwp) :: Chol = .false.
    ! CHOI
    Logical(kind=iwp) :: Choi = .false.
    ! WTHR      thresholds for writing large components in the
    !           first-order perturbed wave function, 3 values that
    !           are for denominator, coefficient, and energy
    Real(kind=wp) :: DnmThr = 0.3_wp,CmpThr = 0.025_wp,CntThr = 0.005_wp
    ! FOCK      string representing the type of Fock matrix
    Character(len=8)  :: FockType = 'STANDARD'
    ! HZER      string representing the type of 0-order hamiltonian
    Character(len=8)  :: Hzero = 'STANDARD'
    ! G1SE      include secondary/inactive elements of the exchange
    !           matrix in the g1 modification to the fock matrix
    Logical(kind=iwp) :: G1SecIn = .false.
    ! RHSD      use the RHS-ondemand algorithm for the calculation of the right-hand side
    Logical(kind=iwp) :: RHSD = .false.
    ! CUMU
    Logical(kind=iwp) :: doCumulant = .false.

  end type ! end of type InputTable

  ! Define the Input as an InputTable structure
  type(InputTable),allocatable :: Input

  save

contains

  subroutine readin_CASPT2(LuIn,nSym)
    !SVC read and store the input as independent as possible. Any sanity
    ! checks not required for reading in the input should be postponed till
    ! the proc_inp call (processing of input). The only variable needed here
    ! is nSym, as some input lines assume knowledge of the number of irreps.

    use definitions, only: iwp,u6

    implicit none

    Integer(kind=iwp),intent(in) :: LuIn,nSym

    Character(len=128) :: Line
    Character(len=:),allocatable :: dLine
    Character(len=4) :: Command,Word

    Integer(kind=iwp) :: i,j,iSym,nStates
    Integer(kind=iwp) :: iSplit,iError

    Logical(kind=iwp),external :: next_non_comment

#ifdef _ENABLE_CHEMPS2_DMRG_
    Logical(kind=iwp) :: dochemps2 = .false.
#endif

    rewind (LuIn)
    call RdNLst(LuIn,'CASPT2')

    ! beginning of reading loop
    do

      if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
      Command = Line(1:4)
      call Upcase(Command)

      !IFG Note that when multiple values are required, ExtendLine may
      ! be called (0 or more times) until the READ statement gives no error
      ! this allows the input to be split in lines more or less arbitrarily,
      ! as if the values were read directly from the file.
      select case (Command)

      case ('TITL')
        if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
        read (Line,'(A128)') Input%Title

        ! File with the reference CAS/RAS wavefunction
      case ('FILE')
        if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
        ! Not using list-directed input (*), because then the slash means end of input
        read (Line,'(A)',IOStat=iError) Input%file
        if (iError /= 0) call IOError(Line)

        ! Root selection
      case ('MULT')
        Input%MULT = .true.
        if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
        read (Line,*) Word
        call Upcase(Word)
        if (Word == 'ALL') then
          nStates = 0
          Input%AllMult = .true.
        else
          read (Line,*,IOStat=iError) nStates
          if (iError /= 0) call IOError(Line)
          if (nStates <= 0) call StatesError(Line)
        end if
        ! TODO: use mma_allocate if possible
        allocate (Input%MultGroup%State(nStates))
        Input%nMultState = nStates
        iSplit = scan(Line,' ')
        allocate (Character(len=len(Line)) :: dLine)
        dLine = Line(iSplit:)
        iError = -1
        do while (iError < 0)
          read (dLine,*,IOStat=iError) (Input%MultGroup%State(i),i=1,nStates)
          if (iError > 0) call IOError(Line)
          if (iError < 0) then
            if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
            call ExtendLine(dLine,Line)
          end if
        end do
        deallocate (dLine)

      case ('XMUL')
        Input%XMUL = .true.
        if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
        read (Line,*) Word
        call Upcase(Word)
        if (Word == 'ALL') then
          nStates = 0
          Input%AllXMult = .true.
        else
          read (Line,*,IOStat=iError) nStates
          if (iError /= 0) call IOError(Line)
          if (nStates <= 0) call StatesError(Line)
        end if
        allocate (Input%XMulGroup%State(nStates))
        Input%nXMulState = nStates
        iSplit = scan(Line,' ')
        allocate (Character(len=len(Line)) :: dLine)
        dLine = Line(iSplit:)
        iError = -1
        do while (iError < 0)
          read (dLine,*,IOStat=iError) (Input%XMulGroup%State(i),i=1,nStates)
          if (iError > 0) call IOError(Line)
          if (iError < 0) then
            if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
            call ExtendLine(dLine,Line)
          end if
        end do
        deallocate (dLine)

      case ('RMUL')
        Input%RMUL = .true.
        if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
        read (Line,*) Word
        call Upcase(Word)
        if (Word == 'ALL') then
          nStates = 0
          Input%AllRMult = .true.
        else
          read (Line,*,IOStat=iError) nStates
          if (iError /= 0) call IOError(Line)
          if (nStates <= 0) call StatesError(Line)
        end if
        allocate (Input%RMulGroup%State(nStates))
        Input%nRMulState = nStates
        iSplit = scan(Line,' ')
        allocate (Character(len=len(Line)) :: dLine)
        dLine = Line(iSplit:)
        iError = -1
        do while (iError < 0)
          read (dLine,*,IOStat=iError) (Input%RMulGroup%State(i),i=1,nStates)
          if (iError > 0) call IOError(Line)
          if (iError < 0) then
            if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
            call ExtendLine(dLine,Line)
          end if
        end do
        deallocate (dLine)

      case ('DWMS')
        Input%DWMS = .true.
        if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
        read (Line,*,IOStat=iError) Input%ZETA
        if (iError /= 0) call IOError(Line)

      case ('DWTY')
        if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
        read (Line,*,IOStat=iError) Input%DWType
        if (iError /= 0) call IOError(Line)

      case ('EFOC')
        Input%EFOC = .true.

      case ('LROO')
        Input%LROO = .true.
        if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
        read (Line,*,IOStat=iError) Input%SingleRoot
        if (iError /= 0) call IOError(Line)

      case ('RLXR')
        if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
        read (Line,*,IOStat=iError) Input%RlxRoot
        if (iError /= 0) call IOError(Line)

        ! freeze-deleted control

      case ('FROZ')
        Input%FROZ = .true.
        allocate (Input%nFro(nSYM))
        if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
        allocate (Character(len=len(Line)) :: dLine)
        dLine = Line
        iError = -1
        do while (iError < 0)
          read (dLine,*,IOStat=iError) (Input%nFro(iSym),iSym=1,nSym)
          if (iError > 0) call IOError(Line)
          if (iError < 0) then
            if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
            call ExtendLine(dLine,Line)
          end if
        end do
        deallocate (dLine)

      case ('DELE')
        Input%DELE = .true.
        allocate (Input%nDel(nSYM))
        if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
        allocate (Character(len=len(Line)) :: dLine)
        dLine = Line
        iError = -1
        do while (iError < 0)
          read (dLine,*,IOStat=iError) (Input%nDel(iSym),iSym=1,nSym)
          if (iError > 0) call IOError(Line)
          if (iError < 0) then
            if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
            call ExtendLine(dLine,Line)
          end if
        end do
        deallocate (dLine)

        ! equation solver control

      case ('MAXI')
        if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
        read (Line,*,IOStat=iError) Input%maxIter
        if (iError /= 0) call IOError(Line)

      case ('CONV')
        if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
        read (Line,*,IOStat=iError) Input%ThrConv
        if (iError /= 0) call IOError(Line)

      case ('THRE')
        Input%THRE = .true.
        if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
        allocate (Character(len=len(Line)) :: dLine)
        dLine = Line
        iError = -1
        do while (iError < 0)
          read (dLine,*,IOStat=iError) Input%ThrsHN,Input%ThrsHS
          if (iError > 0) call IOError(Line)
          if (iError < 0) then
            if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
            call ExtendLine(dLine,Line)
          end if
        end do
        deallocate (dLine)

      case ('SHIF')
        if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
        read (Line,*,IOStat=iError) Input%Shift
        if (iError /= 0) call IOError(Line)

      case ('IMAG')
        if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
        read (Line,*,IOStat=iError) Input%ShiftI
        if (iError /= 0) call IOError(Line)

        ! environment

      case ('RFPE')
        Input%RFpert = .true.

      case ('OFEM')
        Input%OFEmbedding = .true.

        ! print controls

      case ('PRWF')
        if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
        read (Line,*,IOStat=iError) Input%PrWF
        if (iError /= 0) call IOError(Line)

      case ('OUTP')
        if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
        call StdFmt(Line,Input%OutFormat)

      case ('NOOR')
        Input%PrOrb = .false.

      case ('PRSD')
        Input%PRSD = .true.

      case ('WTHR')
        if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
        allocate (Character(len=len(Line)) :: dLine)
        dLine = Line
        iError = -1
        do while (iError < 0)
          read (dLine,*,IOStat=iError) Input%DNMTHR,Input%CMPTHR,Input%CNTTHR
          if (iError > 0) call IOError(Line)
          if (iError < 0) then
            if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
            call ExtendLine(dLine,Line)
          end if
        end do
        deallocate (dLine)

        ! properties

      case ('DENS')
        Input%DENS = .true.

      case ('PROP')
        Input%Properties = .true.

      case ('NOPR')
        Input%Properties = .false.

        ! fock matrix, 0-order hamiltonian

      case ('TRAN')
        Input%ORBIN = 'TRANSFOR'

      case ('FOCK')
        if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
        call StdFmt(Line,Input%FockType)

      case ('HZER')
        if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
        call StdFmt(Line,Input%HZero)

      case ('G1SE')
        Input%G1SECIN = .true.

      case ('IPEA')
        Input%IPEA = .true.
        if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
        read (Line,*,IOStat=iError) Input%BSHIFT
        if (iError /= 0) call IOError(Line)

        ! cholesky

      case ('CHOL')
        Input%Chol = .true.
        call Cho_caspt2_rdInp(.true.,LuIn)

      case ('CHOI')
        Input%ChoI = .true.
        call Cho_caspt2_rdInp(.false.,LuIn)

        ! freeze-delete approximation schemes

      case ('AFRE')
        Input%aFreeze = .true.
        Input%modify_correlating_MOs = .true.
        if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
        allocate (Character(len=len(Line)) :: dLine)
        dLine = Line
        iError = -1
        do while (iError < 0)
          read (dLine,*,IOStat=iError) Input%lnFro,Input%ThrFr,Input%ThrDe
          if (iError > 0) call IOError(Line)
          if (iError < 0) then
            if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
            call ExtendLine(dLine,Line)
          end if
        end do
        deallocate (dLine)
        allocate (Input%NamFro(Input%lnFro))
        if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
        call Upcase(Line)
        allocate (Character(len=len(Line)) :: dLine)
        dLine = Line
        iError = -1
        do while (iError < 0)
          read (dLine,*,IOStat=iError) (Input%NamFro(i),i=1,Input%lnFro)
          if (iError > 0) call IOError(Line)
          if (iError < 0) then
            if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
            call Upcase(Line)
            call ExtendLine(dLine,Line)
          end if
        end do
        deallocate (dLine)

      case ('LOVC')
        Input%LovCASPT2 = .true.
        Input%modify_correlating_MOs = .true.
        if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
        read (Line,*,IOStat=iError) Input%thr_atm
        if (iError /= 0) call IOError(Line)

      case ('FNOC')
        Input%FnoCASPT2 = .true.
        Input%modify_correlating_MOs = .true.
        if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
        read (Line,*,IOStat=iError) Input%vFrac
        if (iError /= 0) call IOError(Line)

      case ('DOMP')
        Input%doMP2 = .true.

      case ('DOEN')
        Input%doEnv = .true.

      case ('VIRA')
        Input%VIRA = .true.

      case ('GHOS')
        Input%GhostDelete = .true.
        Input%modify_correlating_MOs = .true.
        if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
        read (Line,*,IOStat=iError) Input%ThrGD
        if (iError /= 0) call IOError(Line)

      case ('NOMU')
        Input%NoMult = .true.

      case ('ONLY')
        if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
        read (Line,*,IOStat=iError) Input%OnlyRoot
        if (iError /= 0) call IOError(Line)

      case ('NOMI')
        Input%NoMix = .true.

      case ('RHSD')
        Input%RHSD = .true.

#ifdef _ENABLE_BLOCK_DMRG_
      case ('CUMU')
        Input%doCumulant = .true.
#elif _ENABLE_CHEMPS2_DMRG_
      case ('CHEM')
        !Quan: Using the same variable doCumulant in Block
        Input%doCumulant = .true.
        dochemps2 = .true.
#endif

      case ('EFFE')
        Input%JMS = .true.
        if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
        read (Line,*,IOStat=iError) nStates
        if (iError /= 0) call IOError(Line)
        allocate (Input%Heff(nStates,nStates))
        Input%Heff = 0.0_wp  ! 8 bytes
        do i = 1,nStates
          if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
          allocate (Character(len=len(Line)) :: dLine)
          dLine = Line
          iError = -1
          do while (iError < 0)
            read (dLine,*,IOStat=iError) (Input%Heff(i,j),j=1,nStates)
            if (iError > 0) call IOError(Line)
            if (iError < 0) then
              if (.not. next_non_comment(LuIn,Line)) call EOFError(Line)
              call ExtendLine(dLine,Line)
            end if
          end do
          deallocate (dLine)
        end do

        ! OBSOLETE KEYWORDS

      case ('GRAD')
        call WarningMessage(2,'Obsolete keyword: '//Command)
        call Quit_OnUserError

      case ('NOTR')
        call WarningMessage(2,'Obsolete keyword: '//Command)
        call Quit_OnUserError

      case ('JACO')
        call WarningMessage(2,'Obsolete keyword: '//Command)
        call Quit_OnUserError

      case ('EXTR')
        call WarningMessage(2,'Obsolete keyword: '//Command)
        call Quit_OnUserError

      case ('QLQR')
        call WarningMessage(2,'Obsolete keyword: '//Command)
        call Quit_OnUserError

      case ('NATU')
        call WarningMessage(2,'Obsolete keyword: '//Command)
        call Quit_OnUserError

      case ('MOLO')
        call WarningMessage(2,'Obsolete keyword: '//Command)
        call Quit_OnUserError

        ! DONE WITH READING INPUT

      case ('END ')
        exit

        ! NO MATCH FOUND, UNKOWN KEYWORD

      case Default
        call WarningMessage(2,'Unrecognized keyword: '//Command)
        call Quit_OnUserError

      end select

    end do ! end of reading loop

#ifdef _ENABLE_CHEMPS2_DMRG_
    ! Check if nState>1
    if ((dochemps2 .EQV. .true.) .and. (nStates > 1)) then
      write (u6,*) 'CHEMPS2> Only State Specific calculation supported'
      call Quit_OnUserError()
    endif
#endif

    ! Normal exit
    return

  end subroutine readin_CASPT2

  subroutine ExtendLine(DynLine,Line)
    implicit none
    Character(len=:),allocatable,intent(InOut) :: DynLine
    Character(len=*),intent(In)                :: Line
    Character(len=len_trim(DynLine))           :: Aux
    Aux = DynLine
    deallocate (DynLine)
    allocate (Character(len=len(Aux) + len(Line) + 1) :: DynLine)
    DynLine = trim(Aux)//' '//Line
  end subroutine ExtendLine

  subroutine IOError(line)
    use definitions,only:u6
    implicit none
    Character(len=*),intent(in) :: line

    call WarningMessage(2,'I/O error when reading line.')
    write (u6,*) 'Last line read from input: ',line
    call Quit_OnUserError
  end subroutine IOError

  subroutine EOFError(line)
    use definitions,only:u6
    implicit none
    Character(len=*),intent(in) :: line

    call WarningMessage(2,'Premature end of input file.')
    write (u6,*) 'Last line read from input: ',line
    call Quit_OnUserError
  end subroutine EOFError

  subroutine StatesError(line)
    use definitions,only:u6
    implicit none
    Character(len=*),intent(in) :: line

    call WarningMessage(2,'Number of states must be > 0.')
    write (u6,*) 'Last line read from input: ',line
    call Quit_OnUserError
  end subroutine StatesError

end module
