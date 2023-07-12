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
subroutine procinp_caspt2
  !SVC: process CASPT2 input based on the data in the input table, and
  ! initialize global common-block variables appropriately.
  use inputData, only: input
  use definitions, only: iwp,wp
  use caspt2_output, only: iPrGlb, terse, cmpThr, cntThr, dnmThr
  use caspt2_global, only: sigma_p_epsilon, sigma_p_exponent, &
                           ipea_shift, imag_shift, real_shift
  use caspt2_gradient, only: do_grad, do_nac, do_csf, iRoot1, iRoot2
  use UnixInfo, only: SuperName
#ifdef _MOLCAS_MPP_
  use Para_Info, only:Is_Real_Par, nProcs
#endif
! NOT TESTED
#if 0
  use OFembed, only:Do_OFemb
#endif

  implicit none

#include "rasdim.fh"
#include "caspt2.fh"
#include "pt2_guga.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "SysDef.fh"
#include "chocaspt2.fh"
#include "caspt2_grad.fh"

  integer(kind=iwp) :: iDummy

  ! Number of non-valence orbitals per symmetry
  integer(kind=iwp) :: nCore(mxSym)
  integer(kind=iwp) :: nDiff,NFI,NSD
  ! Geometry-determining root
  logical(kind=iwp) :: Is_iRlxRoot_Set, do_real, do_imag, do_sigp
  ! Environment
  character(Len=180) :: Env

  integer(kind=iwp) :: I,J,M,N
  integer(kind=iwp) :: iSym
  ! State selection
  integer(kind=iwp) :: iGroup,iOff
  ! Numerical gradients
  logical(kind=iwp) :: DNG, DNG_available
  integer(kind=iwp) :: iDNG
  integer(kind=iwp), external :: isStructure

  ! Hzero and Focktype are merged together into Hzero. We keep the
  ! variable Focktype not to break the input keyword which is documented
  ! in the manual. However, eventually we will have to keep only Hzero
  ! and remove Focktype.
  Hzero = input%Hzero
  if (Hzero .ne. 'STANDARD' .and. Hzero .ne. 'CUSTOM') then
    call WarningMessage(2,'invalid 0th-order Hamiltonian: '//TRIM(Hzero))
    call Quit_OnUserError
  end if

  ! Choose Focktype, reset IPEA shift to 0 for non-standard fock matrices
  Focktype = input%Focktype
  if (Focktype .ne. 'STANDARD') then
    ! if both Hzero and Focktype are not standard, quit
    if (Hzero .ne. 'STANDARD') then
      call WarningMessage(2,'Requested combination of FOCKtype'//' and HZERo not possible.')
      call Quit_OnUserError
    end if
    ! IPEA shift different from zero only for standard Focktype
    if (ipea_shift .ne. 0.0d0) then
      ipea_shift = 0.0d0
      if (IPRGLB .ge. TERSE) then
        call WarningMessage(1,'IPEA shift reset to zero!')
      end if
    end if
  else
    ! user-specified IPEA shift or not?
    if (input%ipea) then
      ipea_shift = input%ipea_shift
    else
      ! Set default IPEA to 0.25 Eh or 0.0
      call getenvf('MOLCAS_NEW_DEFAULTS',Env)
      call upcase(Env)
      if (Env .eq. 'YES') then
        ipea_shift = 0.0d0
      else
        ipea_shift = 0.25d0
      end if
    end if
  end if

  ! copy over to Hzero the content of Focktype, if Hzero is not CUSTOM
  if (Hzero .ne. 'CUSTOM') then
    Hzero = Focktype
  end if

  ! print warnings if deviating from the default
  if (Hzero .ne. 'STANDARD') then
    call warningmessage(1,'User-modified 0th-order Hamiltonian!')
  end if

  ! real/imaginary shifts
  real_shift = Input%real_shift
  imag_shift = Input%imag_shift

  ! sigma-p regularizers
  if (input%sigma_1_epsilon /= 0.0_wp .and. input%sigma_2_epsilon /= 0.0_wp) then
    call WarningMessage(2,'SIG1 and SIG2 keywords are mutually exclusive')
    call Quit_OnUserError()
  end if

  if (input%sigma_1_epsilon > 0.0_wp) then
    sigma_p_epsilon = Input%sigma_1_epsilon
    sigma_p_exponent = 1
  end if

  if (input%sigma_2_epsilon > 0.0_wp) then
    sigma_p_epsilon = Input%sigma_2_epsilon
    sigma_p_exponent = 2
  end if

  do_real = real_shift > 0.0_wp
  do_imag = imag_shift > 0.0_wp
  do_sigp = sigma_p_epsilon > 0.0_wp
  if ((do_real .and. (do_imag .or. do_sigp)) .or. (do_imag .and. do_sigp)) then
    call WarningMessage(2,'More than one intruder-state removal technique active: &
                           &SHIFt/IMAGinary/SIG1/SIG2 are mutually exclusive!')
    call Quit_OnUserError()
  end if

! RHS algorithm selection
#ifdef _MOLCAS_MPP_
#ifdef _GA_
  ! The RHS on-demand algorithm doesn't handle serial calculations
  ! because it's not adapted for use with regular Work arrays, only
  ! global arrays, and needs to be switched off (using rhsall instead)
  RHSDIRECT = (Is_Real_Par() .AND. Input%RHSD)
#else
  ! Without the Global Arrays library, we can't use the RHSALL2
  ! and ADDRHS algorithms in parallel. Here we force the use of
  ! RHS on-demand instead, depending on if the calculation is
  ! really parallel or not.
  RHSDIRECT = Is_Real_Par()
#endif
#else
  RHSDIRECT = .False.
#endif

  ! Cholesky: set defaults if it was not called during input
  if (.NOT. (Input%ChoI .or. Input%Chol)) then
    call Cho_caspt2_rdInp(.True.,iDummy)
  end if

  !---  Initialize
  IDCIEX = 0
  IEOF1M = 0
  do I = 1,64
    IAD1M(I) = -1
  end do
  RFpert = Input%RFPert

  NTIT = 0
  OUTFMT = 'DEFAULT'
  G1SECIN = .FALSE.
  PRORB = .TRUE.
  PRSD = .FALSE.
  NCASES = 13

  JMS = Input%JMS
  NLYROOT = Input%OnlyRoot
  NLYGROUP = 0

  DoCumulant = Input%DoCumulant

!***********************************************************************
!
! Root selection
!
!***********************************************************************
  !SB: create a flat array, putting (X)MS states with increasing
  ! group number. do not allow to run both MS and XMS in the same
  ! calculation, it can lead to catastrophic results.
  ! For MS, put one state per group, for XMS put all states in
  ! a single group.
  ! Example: lets say the input is:
  ! MULTistate = 3 3 4 5
  !
  ! then the arrays that are created will be:
  ! nstate = 3
  ! mstate(nstate): 3 4 5
  ! nGroup = 3
  ! nGroupState(nGroup): 1 1 1
  !
  ! On the other hand, if the input is
  ! XMULtistate = 4 1 2 3 4
  !
  ! then the arrays that are created will be:
  ! nstate = 4
  ! mstate(nstate): 1 2 3 4
  ! nGroup = 1
  ! nGroupState(nGroup): 4
  NSTATE = 0
  MSTATE = 0
  NGROUP = 0
  NGROUPSTATE = 0
  ! This is the case for MS-CASPT2 and DW-CASPT2
  if (Input%MULT) then
    if (Input%XMUL .or. Input%RMUL) then
      call WarningMessage(2,'Keyword MULT cannot be used with neither XMUL nor RMUL.')
      call Quit_OnUserError
    end if
    ! Either the states were specified manually or the keyword "all"
    ! was used, so first we check the keyword all
    if (Input%AllMult) then
      NSTATE = NROOTS
      MSTATE = IROOT
      NGROUP = NSTATE
      NGROUPSTATE(1:NGROUP) = 1
    else
      ! Save the states that need to be computed
      do I = 1,Input%nMultState
        MSTATE(I) = Input%MultGroup%State(I)
        NSTATE = NSTATE + 1
      end do
      NGROUP = Input%nMultState
      NGROUPSTATE(1:NGROUP) = 1
    end if
  end if
  IOFF = NSTATE

  ! This is the case for XMS-CASPT2 and XDW-CASPT2
  if (Input%XMUL) then
    if (Input%MULT .or. Input%RMUL) then
      call WarningMessage(2,'Keyword XMUL cannot be used with neither MULT nor RMUL.')
      call Quit_OnUserError
    end if
    ! This is a XDW-CASPT2 calculation. It is actually more similar to
    ! a MS-CASPT2 one since we need to put one state per group and thus
    ! have as many groups as states. Nevertheless, it makes more sense
    ! from the user point of view to ask for it through XMUL and DWMS
    if (Input%DWMS) then
      if (Input%AllXMult) then
        NSTATE = NROOTS
        MSTATE = IROOT
        NGROUP = NSTATE
        NGROUPSTATE(1:NGROUP) = 1
      else
        do I = 1,Input%nXMulState
          MSTATE(I) = Input%XMulGroup%State(I)
          NSTATE = NSTATE + 1
        end do
        NGROUP = Input%nXMulState
        NGROUPSTATE(1:NGROUP) = 1
      end if
    ! This is a XMS-CASPT2: one group with all the states
    else
      if (Input%AllXMult) then
        NSTATE = NROOTS
        MSTATE = IROOT
        NGROUP = 1
        NGROUPSTATE(1) = NSTATE
      else
        NGROUP = 1
        NGROUPSTATE(NGROUP) = Input%nXMulState
        do I = 1,Input%nXMulState
          MSTATE(I) = Input%XMulGroup%State(I)
          NSTATE = NSTATE + 1
        end do
      end if
    end if
  end if

  ! This is the case for RMS-CASPT2
  if (Input%RMUL) then
    if (Input%MULT .or. Input%XMUL) then
      call WarningMessage(2,'Keyword RMUL cannot be used with neither MULT nor XMUL.')
      call Quit_OnUserError
    end if
    if (Input%AllRMult) then
      NSTATE = NROOTS
      MSTATE = IROOT
      NGROUP = NSTATE
      NGROUPSTATE(1:NGROUP) = 1
    else
      do I = 1,Input%nRMulState
        MSTATE(I) = Input%RMulGroup%State(I)
        NSTATE = NSTATE + 1
      end do
      NGROUP = Input%nRMulState
      NGROUPSTATE(1:NGROUP) = 1
    end if
  end if
  ! After parsing mult, xmult or rmult, check that no two equal states where
  ! given in the input
  if (Input%MULT .OR. Input%XMUL .or. Input%RMUL) then
    do I = 1,NSTATE
      do J = I + 1,NSTATE
        if (MSTATE(I) .EQ. MSTATE(J)) then
          call WarningMessage(2,'The same root cannot be used twice in MULT/XMUL/RMUL blocks.')
          call Quit_OnUserError
        end if
      end do
    end do
  end if
  ! The LROOt keyword specifies a single root to be used. It should not be
  ! used together with either MULT or XMUL keywords.
  if (Input%LROO) then
    if (Input%MULT .OR. Input%XMUL .or. Input%RMUL) then
      call WarningMessage(2,'Keyword LROO cannot be used together with the MULT or XMUL keywords.')
      call Quit_OnUserError
    end if
    NSTATE = 1
    MSTATE(1) = Input%SingleRoot
    NGROUP = 1
    NGROUPSTATE(1) = 1
  end if
  ! if still nothing was selected we should default to compute all the
  ! roots that were part of the rasscf orbital optimization.
  if (NSTATE .EQ. 0) then
    NSTATE = NROOTS
    MSTATE = IROOT
    NGROUP = NSTATE
    NGROUPSTATE(1:NGROUP) = 1
  end if
  ! Find the group number for OnlyRoot
  if (NLYROOT .ne. 0) then
    IOFF = 0
    do IGROUP = 1,NGROUP
      do I = 1,NGROUPSTATE(IGROUP)
        if (IOFF + I .eq. NLYROOT) NLYGROUP = IGROUP
      end do
      IOFF = IOFF + NGROUPSTATE(IGROUP)
    end do
  end if
  ! Finally, some sanity checks.
  if (NSTATE .LE. 0 .OR. NSTATE .GT. MXROOT) then
    call WarningMessage(2,'Number of states is <0 or too large.')
    write (6,'(a,i8)') ' NSTATE = ',NSTATE
    write (6,*) ' Check usage of keywords MULT/XMUL/RMUL.'
    call Quit_OnUserError
  end if
  ! setup root to state translation
  ROOT2STATE = 0
  do I = 1,NSTATE
    ROOT2STATE(MSTATE(I)) = I
  end do
  !
  ! Relax root selection
  !
  iRlxRoot = -1
  call Qpg_iScalar('NumGradRoot',Is_iRlxRoot_Set)
  if (Is_iRlxRoot_Set) then
    call Get_iScalar('NumGradRoot',iRlxRoot)
  end if
  if (Input%RlxRoot .gt. 0) iRlxRoot = Input%RlxRoot
  if (iRlxRoot .eq. -1) iRlxRoot = NSTATE
  if (iRlxRoot .gt. NSTATE) then
    if (IPRGLB .GE. TERSE) then
      call WarningMessage(1,'Too large iRlxRoot.')
      write (6,*) ' Reset to NSTATE=',NSTATE
    end if
    iRlxRoot = NSTATE
  end if

!***********************************************************************
!
! Determine number of Frozen/Deleted orbitals
!
!***********************************************************************
  ! The number of Frozen orbitals is initially read from the reference
  ! wavefunction. Here, we modify that number to be the larger of what
  ! was in the reference and the non-valence orbitals.
  call Get_iArray('Non valence orbitals',nCore,nSym)
  do iSym = 1,nSym
    if (nCore(iSym) .gt. nFro(iSym)) then
      nDiff = nCore(iSym) - nFro(iSym)
      nDiff = min(nDiff,nISh(iSym))
      nFro(iSym) = nFro(iSym) + nDiff
      nISh(iSym) = nISh(iSym) - nDiff
    end if
  end do
  ! if a user specified the number of frozen orbitals explicitely, then
  ! that number overwrites the automatically chosen numbers, but it we
  ! warn the user and report what the computed number was.
  if (Input%FROZ) then
    if (IPRGLB .GE. TERSE) then
      call WarningMessage(1,'User changed nr of frozen orbitals.')
    end if
    do I = 1,NSYM
      NFI = NFRO(I) + NISH(I)
      if (NFI .LT. Input%nFro(I)) then
        call WarningMessage(2,'Too many frozen orbitals!')
        call Quit_OnUserError
      else
        nFro(I) = Input%nFro(I)
      end if
      NISH(I) = NFI - nFro(I)
    end do
  end if
  ! Set user-specified number of deleted orbitals.
  if (Input%DELE) then
    do I = 1,NSYM
      NSD = NSSH(I) + NDEL(I)
      if (NSD .LT. Input%nDel(I)) then
        call WarningMessage(2,'Too many deleted orbitals!')
        call Quit_OnUserError
      else
        NDEL(I) = Input%nDel(I)
      end if
      NSSH(I) = NSD - nDel(I)
    end do
  end if

! NOT TESTED
#if 0
!***********************************************************************
!
! Orbital-free embedding.
!
!***********************************************************************
  if (Input%OFEmbedding) then
    Do_OFemb = .true.
    write (6,*)
    write (6,*) '  --------------------------------------'
    write (6,*) '   Orbital-Free Embedding Calculation   '
    write (6,*) '  --------------------------------------'
    write (6,*)
  end if
#endif
!***********************************************************************
!
! Determine what kind of calculations are needed/requested.
!
!***********************************************************************
  ! flags for enabling/disabling sections
  IFPROP = Input%Properties
  IFDENS = Input%DENS
  IFMIX = .NOT. Input%NoMix
  IFXMS = Input%XMUL
  IFRMS = Input%RMUL
  IFMSCOUP = (Input%MULT.or.IFXMS.or.IFRMS) .AND. (.NOT. Input%NoMult)
  IF (nState == 1 .AND. (NLYROOT .ne. mState(1))) IFMSCOUP = .FALSE.
  IFDW = Input%DWMS
  ! Set type and exponent for DWMS
  DWType = Input%DWType
  if (IFDW) then
    if (DWType <= 0 .or. DWType > 3) then
      call WarningMessage(2,'DWTYpe should be either 1, 2 or 3.')
      call Quit_OnUserError
    end if
    zeta = Input%zeta
  end if

  ! Choice? of preprocessing route
  ORBIN = 'TRANSFOR'

  !---  Expectation values of H0 needed or not? Spectral resolution?
  BMATRIX = 'YES     '
  BTRANS = 'YES     '
  BSPECT = 'YES     '
  !---  Overlap matrix needed? Linear dependence removal needed?
  SMATRIX = 'YES     '
  SDECOM = 'YES     '
  ! Presently, we always need the ON transformation/Lindep removal.
  ! Consistency of these demands:
  if (BSPECT .NE. 'NO      ') BTRANS = 'YES     '
  if (BTRANS .NE. 'NO      ') BMATRIX = 'YES     '
  if (BTRANS .NE. 'NO      ') SDECOM = 'YES     '
  if (SDECOM .NE. 'NO      ') SMATRIX = 'YES     '

!***********************************************************************
!
! Thresholds
!
!***********************************************************************
  if (Input%THRE) then
    if (IPRGLB .GE. TERSE) then
      call WarningMessage(1,'User modified linear dependency thresholds!')
    end if
  end if
  THRSHN = Input%THRSHN
  THRSHS = Input%THRSHS
  THRCONV = Input%THRCONV
  CITHR = Input%PrWF
  THRENE = 5.0d+01
  THROCC = 5.0d-04
  MAXIT = Input%MaxIter
  DNMTHR = Input%DNMTHR
  CMPTHR = Input%CMPTHR
  CNTTHR = Input%CNTTHR

  !---  Create the symmetry multiplication table
  MUL(1,1) = 1
  M = 1
  do N = 1,3
    do I = 1,M
      do J = 1,M
        MUL(I + M,J) = M + MUL(I,J)
        MUL(I,J + M) = MUL(I + M,J)
        MUL(I + M,J + M) = MUL(I,J)
      end do
    end do
    M = 2*M
  end do
!***********************************************************************
!
! Gradients
!
!***********************************************************************

  ! at the moment the calculation of analytic gradients has many
  ! technical restrictions and we need to make sure that we do not
  ! run into a unsupported combination of keywords, these are:
  ! 1. no ipea shift
  ! 2. no symmetry
  ! 3. all CASSCF roots included in QD-CASPT2 gradients
  ! 4. QD-CASPT2 gradients only with CD/DF

  ! CASPT2 is invariant with respect to rotations in active space?
  INVAR = .true.
  call put_iScalar('mp2prpt',0)
  if (ipea_shift /= 0.0_wp) INVAR = .false.

  ! check if numerical gradients were requested in GATEWAY
  call qpg_iScalar('DNG', DNG_available)
  if (DNG_available) then
    call get_iScalar('DNG', iDNG)
    DNG = iDNG .eq. 1
  else
    DNG = .false.
  end if

  ! if the CASPT2 module was called by the NUMERICAL_GRADIENT one,
  ! we make sure to disable analytic gradients
  if (SuperName(1:18) == 'numerical_gradient') then
    call put_iScalar('mp2prpt',0)
    DNG = .true.
    do_grad = .false.
  end if

  ! check first if the user specifically asked for analytic gradients
  if (input%GRDT) then
    do_grad = Input%GRDT

    ! quit if both analytical and numerical gradients were explicitly requested
    if (DNG) then
      call warningMessage(2,'It seems that numerical gradients were requested'// &
                            ' in GATEWAY and analytical gradients in CASPT2.'// &
                            ' Please choose only one of the two!')
      call quit_onUserError
    end if

    ! only allow analytic gradients without symmetry
    if (nSym /= 1) then
      call warningMessage(2,'Analytic gradients only available without symmetry.')
      call quit_onUserError
    end if

    ! no analytic gradients available with IPEA
    if (.not. INVAR) then
      call warningMessage(2,'Analytic gradients not available with IPEA shift.')
      call quit_onUserError
    end if

    ! only allow analytic gradients either with nstate = nroots or with sadref
    if ((nState /= nRoots) .and. (.not. ifsadref)) then
      call warningMessage(2,'Analytic gradients available only if all'// &
                            ' CASSCF roots are included in the CASPT2'// &
                            ' calculation or with the SADRef keyword.')
      call quit_onUserError
    end if

    ! QD-CASPT2 analytic gradients available only with DF or CD
    if (ifMSCoup .and. (.not. ifChol)) then
      call warningMessage(2,'MS-type analytic gradients available only '//  &
                            'with density fitting or Cholesky decomposition.')
      call quit_onUserError
    end if
#ifdef _MOLCAS_MPP_
    ! for the time being no gradients with MPI
    if (nProcs > 1) then
      call warningMessage(2,'Analytic gradients not available'//  &
                            ' in parallel executions.')
      call quit_onUserError
    end if
#endif

  end if

  ! inside LAST_ENERGY we do not need analytic gradients
  if (SuperName(1:11) == 'last_energy') DNG=.true.

  ! check if the calculation is inside a loop and make analytical
  ! gradients default in this case, unless the user specifically
  ! requested numerical gradients in GATEWAY
  if ((isStructure() == 1)) then
    ! if MPI is enabled, analytic gradients only with one process
#ifdef _MOLCAS_MPP_
    if (nProcs == 1) then
#endif
      ! check the hard constraints first
      if ((.not. DNG) .and. (nSym == 1) .and. INVAR) then
        do_grad = .true.

        ! check weaker constraints, if not met, revert to numerical gradients
        if (ifMSCoup .and. (.not. ifChol)) do_grad = .false.
        if ((nState /= nRoots) .and. (.not. ifsadref)) do_grad = .false.
      end if
#ifdef _MOLCAS_MPP_
    end if
#endif
  end if

  ! compute full unrelaxed density for gradients
  if (do_grad) ifDens = .true.

  if (do_grad) then
    call put_iScalar('mp2prpt',2)
    do_nac = input%NAC
    if (input%iNACRoot1 == 0 .and. input%iNACRoot2 == 0) then
      iRoot1 = iRlxRoot
      iRoot2 = iRlxRoot
    else
      iRoot1 = input%iNACRoot1
      iRoot2 = input%iNACRoot2
    end if
  end if

  if (do_nac) then
    do_csf = input%CSF
  end if

  IFSADREF = input%SADREF
  IFDORTHO = input%DORTHO

  !! Whether the Fock matrix (eigenvalues) is constructed with
  !! the state-averaged density matrix or not.
  !! The name of the variable is like state-specific DM,
  !! but not necessarily state-specific. It is a matter of the
  !! structure of WORK(LDWGT) array or matrix.
  !! WORK(LDWGT) is a matrix form for SS- and MS-CASPT2 with
  !! state-specific DM, XDW-CASPT2, and RMS-CASPT2, while it is an
  !! array for SS- and MS-CASPT2 with state-averaged DM (with SADREF
  !! option) and XMS-CASPT2.
  if (IFSADREF .or. nRoots == 1 .or. (IFXMS .and. .not.IFDW)) then
    IFSSDM = .false.
  else
    IFSSDM = .true.
  end if

end subroutine procinp_caspt2
