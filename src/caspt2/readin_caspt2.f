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
      Module InputData
      Implicit None
!SVC: this module contains a data structure to keep all input variables.
! Note that I use the standard 'allocate' here against the appropriate
! Molcas practices. The reason is that these are (i) very small, and
! (ii) there is need for allocating complex things such as derived
! types, which are not supported with stdalloc. Hence, the infraction.

#include "compiler_features.h"

      Type States
        Integer, Allocatable :: State(:)
      End Type

      Type InputTable
!       TITL      one line with a descriptive name
        Character(Len=128) :: Title = ' '
!       FILE      file to read CAS/RAS reference from
        Character(Len=128) :: File = 'JOBIPH'
!       MULT      the number of states, followed by the ID of each state
        Logical :: MULT = .False.
        Integer :: nMultState = 0
        Type(States) :: MultGroup
        Logical :: AllMult = .False.
!       XMUL      extended multi-state caspt2
        Logical :: XMUL = .False.
        Integer :: nXMulState = 0
        Type(States) :: XMulGroup
        Logical :: AllXMult = .False.
!       RMUL      rotated multi-state caspt2
        Logical :: RMUL = .False.
        Integer :: nRMulState = 0
        Type(States) :: RMulGroup
        Logical :: AllRMult = .False.
!       DWMS      use dynamical weighting to construct Fock
        Logical :: DWMS = .False.
        Integer :: ZETA = 50
!       EFOC      uses rotated E_0 energies with DWMS
        Logical :: EFOC = .False.
!       LROO      compute only a single root, mutually exclusive
!                 with both MULT or XMUL
        Logical :: LROO = .False.
        Integer :: SingleRoot = 0
!       RLXR      root for which the gradient is computed
        Integer :: RlxRoot = -1

!       IPEA      sets the IP-EA shift
        Logical :: IPEA = .False.
        Real*8 :: BSHIFT = 0.0d0

!       IMAG      size of extra 'imaginary' denominator shift
        Real*8 :: ShiftI = 0.0d0
!       SHIF      size of extra denominator shift
        Real*8 :: Shift = 0.0d0

!       several freeze-delete schemes, each of these should active
!       the general flag below, to indicate additional conversion is
!       needed on the input orbitals
        Logical :: modify_correlating_MOs = .False.
!       AFRE      freeze orbitals that do not have sufficient density
!                 on specified 'active' atoms
        Logical :: aFreeze = .False.
        Integer :: lnFro = 0
        Real*8 :: ThrFr = 0.0d0, ThrDe = 0.0d0
        Character(Len=4), Allocatable :: NamFro(:)
!       LOVC      freeze orbitals that are not localized no the active
!                 site
        Logical :: LovCASPT2 = .False.
        Real*8 :: Thr_Atm = 0.0d0
!       FNOC      delete a fraction of virtual orbitals
        Logical :: FnoCASPT2 = .False.
        Real*8 :: VFrac = 0.0d0
!       DOMP
        Logical :: DoMP2 = .False.
!       DOEN
        Logical :: DoEnv = .False.
!       VIRA
        Logical :: VIRA = .False.
!       GHOS      excludes ghost orbitals from the PT2 treatment
        Logical :: GhostDelete = .False.
        Real*8 :: ThrGD = 0.0d0

!       FROZ      number of frozen orbitals in each irrep
        Logical :: FROZ = .False.
        Integer, Allocatable :: nFro(:)
!       DELE      number of deleted orbitals in each irrep
        Logical :: DELE = .False.
        Integer, Allocatable :: nDel(:)
!       DENS      computes full density matrix from the 1st-order
!                 wavefunction
        Logical :: DENS = .False.
!       RFPE      make a perturbative reaction field calculation
        Logical :: RFPert = .False.
!       THRE      thresholds for removal of:
!         ThrsHN    zero-norm components in the first-order perturbed
!                   wave function
!         ThrsHS    linear dependencies between components of the first-
!                   order perturbed wave function
        Logical :: THRE = .False.
        Real*8 :: ThrsHN = 1.0d-10, ThrsHS = 1.0d-8
!       MAXI      maximum number of iterations for solving a system of
!                 linear equations, default 20. A 0 indicates: use of
!                 the diagonal zeroth order hamiltonian
        Integer :: maxIter = 20
!       Conv      convergence criteria for solving a system of linear
!                 equations
        Real*8 :: ThrConv = 1.0d-6
!       NOMI      do not create an PM-CAS wavefunction file (JobMix)
        Logical :: NoMix = .False.
!       NOMU      do not perform a multistate interaction
        Logical :: noMult = .False.
!       ONLY      in a MS calculation, compute a single root with
!                 couplings to the other roots
        Integer :: OnlyRoot = 0
!       EFFE      read HEff coupling terms from the input and perform
!                 only the multistate part
        Logical :: JMS = .False.
        Real*8, Allocatable :: HEff(:,:)
!       NOOR      do not print orbitals
        Logical :: PrOrb = .True.
!       PROP      compute properties
!       NOPR      do not compute properties
        Logical :: Properties = .False.
!       transformation of reference (input) orbitals
!       NOTR      do not transform to quasi-canonical orbitals,
!                 regardless of the state of the reference orbitals
!       TRAN      transform to quasi-canonical orbitals, regardless
!                 of the state of the reference orbitals
!       the default is to use transformation, unless the PT2 keyword
!       was used in the rasscf program and the fock matrix is standard
        Character(Len=8) :: ORBIN = 'TRANSFOR'
!       OFEM      add orbital-free embedding potential to hamiltonian
        Logical :: OFEmbedding = .False.
!       OUTP      control extent of orbital printing
        Character(Len=8) :: OutFormat = 'DEFAULT '
!       PRWF      print the CI coefficients above this threshold
        Real*8 :: PrWF = 0.05d0
!       PRSD      print the determinant expansion of CSFs
        Logical :: PRSD = .False.
!       NOOR      do not print any orbitals
        Logical :: NoOrb = .False.

!       UNDOCUMENTED KEYWORDS
!       CHOL
        Logical :: CHOL = .False.
!       CHOI
        Logical :: CHOI = .False.
!       WTHR      thresholds for writing large components in the
!                 first-order perturbed wave function, 3 values that
!                 are for denominator, coefficient, and energy
        Real*8 :: DNMTHR = 0.3d0, CMPTHR = 0.025d0, CNTTHR = 0.005d0
!       FOCK      string representing the type of Fock matrix
        Character(Len=8) :: FockType = 'STANDARD'
!       HZER      string representing the type of 0-order hamiltonian
        Character(Len=8) :: HZero = 'STANDARD'
!       G1SE      include secondary/inactive elements of the exchange
!                 matrix in the g1 modification to the fock matrix
        Logical :: G1SecIn = .False.
!       RHSD      use the RHS-ondemand algorithm for the calculation of
!                 the right-hand side
        Logical :: RHSD = .False.
!       CUMU
        Logical :: DoCumulant = .False.

      End Type ! end of type InputTable

      ! Define the Input as an InputTable structure
      Type(InputTable), Allocatable :: Input

      Save

      Contains


      Subroutine Readin_CASPT2(LuIn,nSym)
!SVC Read and store the input as independent as possible. Any sanity
! checks not required for reading in the input should be postponed till
! the proc_inp call (processing of input). The only variable needed here
! is nSym, as some input lines assume knowledge of the number of irreps.

      use definitions, only:u6

      Implicit None

      Integer, intent(in) :: LuIn, nSym

      Character(Len=128) :: Line
      Character(len=:), Allocatable :: dLine
      Character(Len=4) :: Command, Word

      Integer :: i, j, iSym, nStates
      Integer :: iSplit, iError

      logical, external :: next_non_comment

#ifdef _ENABLE_CHEMPS2_DMRG_
      logical :: dochemps2 = .false.
#endif

      Rewind(LuIn)
      Call RdNLst(LuIn,'CASPT2')

      ! TODO: replace this by a while cycle
10    Continue

      if (.not.next_non_comment(LuIn,Line)) call EOFError(Line)
      Command = Line(1:4)
      call UpCase(Command)

!IFG Note that when multiple values are required, ExtendLine may
! be called (0 or more times) until the READ statement gives no error
! this allows the input to be split in lines more or less arbitrarily,
! as if the values were read directly from the file.
      Select Case (Command)

      Case('TITL')
      if (.not.next_non_comment(LuIn,Line)) then
        call EOFError(Line)
      end if
      Read (Line,'(A128)') Input%Title

      ! File with the reference CAS/RAS wavefunction
      Case('FILE')
      if (.not.next_non_comment(LuIn,Line)) then
        call EOFError(Line)
      end if
! Not using list-directed input (*), because then the slash means
! end of input
      read(Line,'(A)',IOStat=iError) Input%file
      if (iError /= 0 ) then
        call IOError(Line)
      end if

      ! Root selection
      Case('MULT')
      Input%MULT = .True.
      if (.not.next_non_comment(LuIn,Line)) then
        call EOFError(Line)
      end if
      Read(Line,*) Word
      Call UpCase(Word)
      ! Either compute all states
      If (Word == 'ALL') Then
        nStates = 0
        Input%AllMult = .True.
      ! or read how many states
      else
        read(Line,*,IOStat=iError) nStates
        if (iError /= 0 ) then
          call IOError(Line)
        end if
        If (nStates <= 0) Then
          call WarningMessage(2,'Number of MULT states must be > 0.')
          write(u6,*)'Last line read from input: ', Line
          call Quit_OnUserError
        End If
      End If
      ! TODO: use mma_allocate if possible
      Allocate(Input%MultGroup%State(nStates))
      Input%nMultState = nStates
      iSplit = SCAN(Line,' ')
      Allocate(Character(Len=Len(Line)) :: dLine)
      dLine = Line(iSplit:)
      iError = -1
      Do While (iError < 0)
        Read(dLine,*,IOStat=iError)
     &    (Input%MultGroup%State(i), i=1,nStates)
        If (iError > 0) then
          call IOError(Line)
        end if
        If (iError < 0) Then
          if (.not.next_non_comment(LuIn,Line)) then
            call EOFError(Line)
          end if
          Call ExtendLine(dLine,Line)
        End If
      End Do
      Deallocate(dLine)

      Case('XMUL')
      Input%XMUL = .True.
      if (.not.next_non_comment(LuIn,Line)) then
        call EOFError(Line)
      end if
      Read(Line,*) Word
      Call UpCase(Word)
      If (Word == 'ALL') Then
        nStates = 0
        Input%AllXMult = .True.
      Else
        Read(Line,*,IOStat=iError) nStates
        if (iError /= 0 ) then
          call IOError(Line)
        end if
        If (nStates <= 0) Then
          call WarningMessage(2,'Number of XMUL states must be > 0.')
          write(u6,*)'Last line read from input: ', Line
          call Quit_OnUserError
        End If
      End If
      Allocate(Input%XMulGroup%State(nStates))
      Input%nXMulState = nStates
      iSplit = SCAN(Line,' ')
      Allocate(Character(Len=Len(Line)) :: dLine)
      dLine = Line(iSplit:)
      iError = -1
      Do While (iError < 0)
        Read(dLine,*,IOStat=iError)
     &      (Input%XMulGroup%State(i), i=1,nStates)
        If (iError > 0) then
          call IOError(Line)
        end if
        If (iError < 0) Then
          if (.not.next_non_comment(LuIn,Line)) then
            call EOFError(Line)
          end if
          Call ExtendLine(dLine,Line)
        End If
      End Do
      Deallocate(dLine)

      Case('RMUL')
      Input%RMUL = .True.
      if (.not.next_non_comment(LuIn,Line)) then
        call EOFError(Line)
      end if
      Read(Line,*) Word
      Call UpCase(Word)
      If (Word=='ALL') Then
        nStates = 0
        Input%AllRMult = .True.
      Else
        Read(Line,*,IOStat=iError) nStates
        if (iError /= 0 ) then
          call IOError(Line)
        end if
        If (nStates <= 0) Then
          call WarningMessage(2,'Number of RMUL states must be > 0.')
          write(u6,*)'Last line read from input: ', Line
          call Quit_OnUserError
        End If
      End If
      Allocate(Input%RMulGroup%State(nStates))
      Input%nRMulState = nStates
      iSplit = SCAN(Line,' ')
      Allocate(Character(Len=Len(Line)) :: dLine)
      dLine = Line(iSplit:)
      iError = -1
      Do While (iError < 0)
        Read(dLine,*,IOStat=iError)
     &      (Input%RMulGroup%State(i), i=1,nStates)
        If (iError > 0) then
          call IOError(Line)
        end if
        If (iError < 0) Then
          if (.not.next_non_comment(LuIn,Line)) then
            call EOFError(Line)
          end if
          Call ExtendLine(dLine,Line)
        End If
      End Do
      Deallocate(dLine)


      Case('DWMS')
      Input % DWMS = .True.
      if (.not.next_non_comment(LuIn,Line)) then
        call EOFError(Line)
      end if
      Read(Line,*,IOStat=iError) Input % ZETA
      if (iError /= 0 ) then
        call IOError(Line)
      end if

      Case('EFOC')
      Input % EFOC = .True.

      Case('LROO')
      Input % LROO = .True.
      if (.not.next_non_comment(LuIn,Line)) then
        call EOFError(Line)
      end if
      Read(Line,*,IOStat=iError) Input % SingleRoot
      if (iError /= 0 ) then
        call IOError(Line)
      end if

      Case('RLXR')
      if (.not.next_non_comment(LuIn,Line)) then
        call EOFError(Line)
      end if
      Read(Line,*,IOStat=iError) Input % RlxRoot
      if (iError /= 0 ) then
        call IOError(Line)
      end if

      ! freeze-deleted control

      Case('FROZ')
      Input % FROZ = .True.
      Allocate(Input % nFro(nSYM))
      if (.not.next_non_comment(LuIn,Line)) then
        call EOFError(Line)
      end if
      Allocate(Character(Len=Len(Line)) :: dLine)
      dLine = Line
      iError = -1
      Do While (iError < 0)
        Read(dLine,*,IOStat=iError) (Input % nFro(iSym), iSym=1,nSym)
        If (iError > 0) then
          call IOError(Line)
        end if
        If (iError < 0) Then
          if (.not.next_non_comment(LuIn,Line)) then
            call EOFError(Line)
          end if
          Call ExtendLine(dLine,Line)
        End If
      End Do
      Deallocate(dLine)

      Case('DELE')
      Input % DELE = .True.
      Allocate(Input % nDel(nSYM))
      if (.not.next_non_comment(LuIn,Line)) then
        call EOFError(Line)
      end if
      Allocate(Character(Len=Len(Line)) :: dLine)
      dLine = Line
      iError = -1
      Do While (iError < 0)
        Read(dLine,*,IOStat=iError) (Input % nDel(iSym), iSym=1,nSym)
        If (iError > 0) then
          call IOError(Line)
        end if
        If (iError < 0) Then
          if (.not.next_non_comment(LuIn,Line)) then
            call EOFError(Line)
          end if
          Call ExtendLine(dLine,Line)
        End If
      End Do
      Deallocate(dLine)

      ! equation solver control

      Case('MAXI')
      if (.not.next_non_comment(LuIn,Line)) then
        call EOFError(Line)
      end if
      Read(Line,*,IOStat=iError) Input % maxIter
      if (iError /= 0 ) then
        call IOError(Line)
      end if

      Case('CONV')
      if (.not.next_non_comment(LuIn,Line)) then
        call EOFError(Line)
      end if
      Read(Line,*,IOStat=iError) Input % ThrConv
      if (iError /= 0 ) then
        call IOError(Line)
      end if

      Case('THRE')
      Input % THRE = .True.
      if (.not.next_non_comment(LuIn,Line)) then
        call EOFError(Line)
      end if
      Allocate(Character(Len=Len(Line)) :: dLine)
      dLine = Line
      iError = -1
      Do While (iError < 0)
        Read(dLine,*,IOStat=iError) Input % ThrsHN, Input % ThrsHS
        If (iError > 0) then
          call IOError(Line)
        end if
        If (iError < 0) Then
          if (.not.next_non_comment(LuIn,Line)) then
            call EOFError(Line)
          end if
          Call ExtendLine(dLine,Line)
        End If
      End Do
      Deallocate(dLine)

      Case('SHIF')
      if (.not.next_non_comment(LuIn,Line)) then
        call EOFError(Line)
      end if
      Read(Line,*,IOStat=iError) Input % Shift
      if (iError /= 0 ) then
        call IOError(Line)
      end if

      Case('IMAG')
      if (.not.next_non_comment(LuIn,Line)) then
        call EOFError(Line)
      end if
      Read(Line,*,IOStat=iError) Input % ShiftI
      if (iError /= 0 ) then
        call IOError(Line)
      end if

      ! environment

      Case('RFPE')
      Input % RFpert=.true.

      Case('OFEM')
      Input % OFEmbedding = .True.
c      call Quit_OnInstError

      ! print controls

      Case('PRWF')
      if (.not.next_non_comment(LuIn,Line)) then
        call EOFError(Line)
      end if
      Read(Line,*,IOStat=iError) Input % PrWF
      if (iError /= 0 ) then
        call IOError(Line)
      end if

      Case('OUTP')
      if (.not.next_non_comment(LuIn,Line)) then
        call EOFError(Line)
      end if
      Call StdFmt(Line,Input % OutFormat)

      Case('NOOR')
      Input % PrOrb = .False.

      Case('PRSD')
      Input % PRSD = .True.

      Case('WTHR')
      if (.not.next_non_comment(LuIn,Line)) then
        call EOFError(Line)
      end if
      Allocate(Character(Len=Len(Line)) :: dLine)
      dLine = Line
      iError = -1
      Do While (iError < 0)
        Read(dLine,*,IOStat=iError)
     &    Input % DNMTHR, Input % CMPTHR, Input % CNTTHR
        If (iError > 0) then
          call IOError(Line)
        end if
        If (iError < 0) Then
          if (.not.next_non_comment(LuIn,Line)) then
            call EOFError(Line)
          end if
          Call ExtendLine(dLine,Line)
        End If
      End Do
      Deallocate(dLine)


      ! properties

      Case('DENS')
      Input % DENS = .True.

      Case('PROP')
      Input % Properties = .True.

      Case('NOPR')
      Input % Properties = .False.

      ! fock matrix, 0-order hamiltonian

      Case('TRAN')
      Input % ORBIN = 'TRANSFOR'

      Case('FOCK')
      if (.not.next_non_comment(LuIn,Line)) then
        call EOFError(Line)
      end if
      Call StdFmt(Line,Input % FockType)

      Case('HZER')
      if (.not.next_non_comment(LuIn,Line)) then
        call EOFError(Line)
      end if
      Call StdFmt(Line,Input % HZero)

      Case('G1SE')
      Input % G1SECIN = .True.

      Case('IPEA')
      Input % IPEA = .True.
      if (.not.next_non_comment(LuIn,Line)) then
        call EOFError(Line)
      end if
      Read(Line,*,IOStat=iError) Input % BSHIFT
      if (iError /= 0 ) then
        call IOError(Line)
      end if

      ! cholesky

      Case('CHOL')
      Input % Chol = .True.
      Call Cho_caspt2_rdInp(.True.,LuIn)

      Case('CHOI')
      Input % ChoI = .True.
      Call Cho_caspt2_rdInp(.False.,LuIn)

      ! freeze-delete approximation schemes

      Case('AFRE')
      Input % aFreeze = .True.
      Input % modify_correlating_MOs = .True.
      if (.not.next_non_comment(LuIn,Line)) then
        call EOFError(Line)
      end if
      Allocate(Character(Len=Len(Line)) :: dLine)
      dLine = Line
      iError = -1
      Do While (iError < 0)
        Read(dLine,*,IOStat=iError)
     &    Input % lnFro, Input % ThrFr, Input % ThrDe
        If (iError > 0) then
          call IOError(Line)
        end if
        If (iError < 0) Then
          if (.not.next_non_comment(LuIn,Line)) then
            call EOFError(Line)
          end if
          Call ExtendLine(dLine,Line)
        End If
      End Do
      Deallocate(dLine)
      Allocate(Input % NamFro(Input % lnFro))
      if (.not.next_non_comment(LuIn,Line)) then
        call EOFError(Line)
      end if
      Call UpCase(Line)
      Allocate(Character(Len=Len(Line)) :: dLine)
      dLine = Line
      iError = -1
      Do While (iError < 0)
        Read(dLine,*,IOStat=iError)
     &    (Input%NamFro(i), i=1,Input%lnFro)
        If (iError > 0) then
          call IOError(Line)
        end if
        If (iError < 0) Then
          if (.not.next_non_comment(LuIn,Line)) then
            call EOFError(Line)
          end if
          Call UpCase(Line)
          Call ExtendLine(dLine,Line)
        End If
      End Do
      Deallocate(dLine)

      Case('LOVC')
      Input % LovCASPT2 = .True.
      Input % modify_correlating_MOs = .True.
      if (.not.next_non_comment(LuIn,Line)) then
        call EOFError(Line)
      end if
      Read(Line,*,IOStat=iError) Input % thr_atm
      if (iError /= 0 ) then
        call IOError(Line)
      end if

      Case('FNOC')
      Input % FnoCASPT2 = .True.
      Input % modify_correlating_MOs = .True.
      if (.not.next_non_comment(LuIn,Line)) then
        call EOFError(Line)
      end if
      Read(Line,*,IOStat=iError) Input % vFrac
      if (iError /= 0 ) then
        call IOError(Line)
      end if

      Case('DOMP')
      Input % DoMP2 = .True.

      Case('DOEN')
      Input % DoEnv = .True.

      Case('VIRA')
      Input % VIRA = .True.

      Case('GHOS')
      Input % GhostDelete = .True.
      Input % modify_correlating_MOs = .True.
      if (.not.next_non_comment(LuIn,Line)) then
        call EOFError(Line)
      end if
      Read(Line,*,IOStat=iError) Input % ThrGD
      if (iError /= 0 ) then
        call IOError(Line)
      end if

      Case('NOMU')
      Input % NoMult = .True.

      Case('ONLY')
      if (.not.next_non_comment(LuIn,Line)) then
        call EOFError(Line)
      end if
      Read(Line,*,IOStat=iError) Input % OnlyRoot
      if (iError /= 0 ) then
        call IOError(Line)
      end if

      Case('NOMI')
      Input % NoMix = .True.

      Case('RHSD')
      Input % RHSD = .True.

#ifdef _ENABLE_BLOCK_DMRG_
      Case('CUMU')
      Input % DoCumulant = .True.
#elif _ENABLE_CHEMPS2_DMRG_
      Case('CHEM')
!Quan: Using the same variable DoCumulant in Block
      Input % DoCumulant = .True.
      dochemps2 = .True.
!      if (nStates > 1) then
!       write(6,*) 'CHEMPS2> Only State Specific calculation supported'
!Quan: FIXME: nStates not defined
!       Call Quit_OnUserError
!      endif
#endif

      Case('EFFE')
      Input % JMS = .True.
      if (.not.next_non_comment(LuIn,Line)) then
        call EOFError(Line)
      end if
      Read(Line,*,IOStat=iError) nStates
      if (iError /= 0 ) then
        call IOError(Line)
      end if
      Allocate(Input % HEff(nStates,nStates))
      Input % HEff = 0.0d0
      Do i=1,nStates
        if (.not.next_non_comment(LuIn,Line)) then
          call EOFError(Line)
        end if
        Allocate(Character(Len=Len(Line)) :: dLine)
        dLine = Line
        iError = -1
        Do While (iError < 0)
          Read(dLine,*,IOStat=iError)
     &      (Input % HEff(i,j),j=1,nStates)
          If (iError > 0) then
            call IOError(Line)
          end if
          If (iError < 0) Then
            if (.not.next_non_comment(LuIn,Line)) then
              call EOFError(Line)
            end if
            Call ExtendLine(dLine,Line)
          End If
        End Do
        Deallocate(dLine)
      End Do


      ! OBSOLETE KEYWORDS

      Case('GRAD')
      Call WarningMessage(2,'Obsolete keyword: '//Command)
      Call Quit_OnUserError

      Case('NOTR')
      Call WarningMessage(2,'Obsolete keyword: '//Command)
      Call Quit_OnUserError

      Case('JACO')
      Call WarningMessage(2,'Obsolete keyword: '//Command)
      Call Quit_OnUserError

      Case('EXTR')
      Call WarningMessage(2,'Obsolete keyword: '//Command)
      Call Quit_OnUserError

      Case('QLQR')
      Call WarningMessage(2,'Obsolete keyword: '//Command)
      Call Quit_OnUserError

      Case('NATU')
      Call WarningMessage(2,'Obsolete keyword: '//Command)
      Call Quit_OnUserError

      Case('MOLO')
      Call WarningMessage(2,'Obsolete keyword: '//Command)
      Call Quit_OnUserError

      ! DONE WITH READING INPUT

      Case('END ')
      GoTo 9000

      ! NO MATCH FOUND, UNKOWN KEYWORD

      Case Default
      Call WarningMessage(2,'Unrecognized keyword: '//Command)
      Call Quit_OnUserError

      End Select
      GoTo 10

9000  CONTINUE

#ifdef _ENABLE_CHEMPS2_DMRG_
! Check if nState>1
      if ((dochemps2.EQV..True.) .and. (nStates > 1)) then
        write(6,*) 'CHEMPS2> Only State Specific calculation supported'
        Call Quit_OnUserError()
      endif
#endif

!---  Normal exit
      Return


      End Subroutine readin_CASPT2


      subroutine ExtendLine(DynLine,Line)
        Implicit None
        Character(Len=:), Allocatable, Intent(InOut) :: DynLine
        Character(Len=*), Intent(In) :: Line
        Character(Len=Len_Trim(DynLine)) :: Aux
        Aux = DynLine
        Deallocate(DynLine)
        Allocate(Character(Len=Len(Aux)+Len(Line)+1) :: DynLine)
        DynLine = Trim(Aux) // ' ' // Line
      end subroutine ExtendLine


      subroutine IOError(line)
        use definitions, only:u6
        implicit none
        character(len=*), intent(in) :: line

        call IOError(Line)
        write(u6,*)'Last line read from input: ', line
        call Quit_OnUserError

      end subroutine IOError


      subroutine EOFError(line)
        use definitions, only:u6
        implicit none
        character(len=*), intent(in) :: line

        call EOFError(Line)
        write(u6,*)'Last line read from input: ', line
        call Quit_OnUserError

      end subroutine EOFError


      End Module
