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

subroutine procinp_caspt2()
!SVC: process CASPT2 input based on the data in the input table, and
! initialize global common-block variables appropriately.

use inputData, only: input
use PrintLevel, only: TERSE
use UnixInfo, only: SuperName
use Molcas, only: MxRoot, MxSym
use caspt2_global, only: cmpThr, cntThr, CompressMPS, ConvInvar, dnmThr, do_csf, do_grad, do_lindep, do_nac, if_equalW, if_invar, &
                         if_invaria, if_SSDM, imag_shift, iParRHS, ipea_shift, iPrGlb, iRoot1, iRoot2, MAXBUF, real_shift, &
                         sigma_p_epsilon, sigma_p_exponent, Weight
use caspt2_module, only: BMatrix, BSpect, BTrans, CIThr, CPT2Method, DMRG, DoCumulant, DWType, FockType, G1SECIN, HZero, IfChol, &
                         IfDens, IfDOrtho, IfDW, IfMix, IFMSCoup, IfProp, IfRMS, IfsadRef, IfXMS, iRlxRoot, iRoot, JMS, MaxIt, &
                         mState, nCases, nDel, nFro, nGroup, nGroupState, nIsh, nLYGroup, nLYRoot, nRas1T, nRas3T, nRoots, nRoots, &
                         nSsh, nState, nSym, OrbIn, OutFmt, PrOrb, PRSD, PT2Method, RFPERT, RHSDirect, Root2State, SDECOM, &
                         SMatrix, ThrConv, ThrEne, ThrOCC, ThrSHN, ThrSHS, Zeta
use SC_NEVPT2, only: Do_FIC, Do_SC, SC_amplitude, SC_prop, SC_thres
#ifdef _DMRG_
use qcmaquis_info, only: qcm_group_names
use qcmaquis_interface_cfg, only: dmrg_file, qcmaquis_param
use qcmaquis_interface, only: qcmaquis_interface_init_checkpoint, qcmaquis_interface_remove_param, qcmaquis_interface_set_param
use PrintLevel, only: VERBOSE
use stdalloc, only: mma_allocate
#endif
#if 0
! NOT TESTED
use OFembed, only: Do_OFemb
#endif
use Constants, only: Zero, Quart
use Definitions, only: wp, iwp, u6, RtoB
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par, nProcs
use Definitions, only: MPIInt
#endif

implicit none
integer(kind=iwp) :: I, iDNG, iDummy, iGroup, iOff, iSym, J, nCore(mxSym), nDiff, NFI, NSD
logical(kind=iwp) :: DNG, DNG_available, do_imag, do_real, do_sigp, Found, Is_iRlxRoot_Set
character(len=180) :: Env
character(len=16) :: mstate1
integer(kind=iwp), external :: isStructure
logical(kind=iwp), external :: RF_On

! nCore: Number of non-valence orbitals per symmetry
! Is_iRlxRoot_Set: Geometry-determining root
! Env: Environment
! mstate1: NAC or not
! iGroup: State selection
! DNG: Numerical gradients

! Hzero and Focktype are merged together into Hzero. We keep the
! variable Focktype not to break the input keyword which is documented
! in the manual. However, eventually we will have to keep only Hzero
! and remove Focktype.
Hzero = input%Hzero
if ((Hzero /= 'STANDARD') .and. (Hzero /= 'CUSTOM') .and. (Hzero /= 'DYALL')) then
  call WarningMessage(2,'invalid 0th-order Hamiltonian: '//trim(Hzero))
  call Quit_OnUserError()
end if

! Choose Focktype, reset IPEA shift to 0 for non-standard fock matrices
Focktype = input%Focktype
if (Focktype /= 'STANDARD') then
  ! if both Hzero and Focktype are not standard, quit
  if ((Hzero /= 'STANDARD') .and. (Hzero /= 'DYALL')) then
    call WarningMessage(2,'Requested combination of FOCKtype'//' and HZERo not possible.')
    call Quit_OnUserError()
  end if
  ! IPEA shift different from zero only for standard Focktype
  if (ipea_shift /= Zero) then
    ipea_shift = Zero
    if (IPRGLB >= TERSE) call WarningMessage(1,'IPEA shift reset to zero!')
  end if
else if (Hzero /= 'DYALL') then
  ! user-specified IPEA shift or not?
  if (input%ipea) then
    ipea_shift = input%ipea_shift
  else
    ! Set default IPEA to 0.25 Eh or 0.0
    call getenvf('MOLCAS_NEW_DEFAULTS',Env)
    call upcase(Env)
    if (Env == 'YES') then
      ipea_shift = Zero
    else
      ipea_shift = Quart
    end if
  end if
end if

! copy over to Hzero the content of Focktype, if Hzero is not CUSTOM or DYALL
if ((Hzero /= 'CUSTOM') .and. (Hzero /= 'DYALL')) Hzero = Focktype

! print warnings if deviating from the default
if ((Hzero /= 'STANDARD') .and. (Hzero /= 'DYALL')) call warningmessage(1,'User-modified 0th-order Hamiltonian!')

! real/imaginary shifts
real_shift = Input%real_shift
imag_shift = Input%imag_shift

! sigma-p regularizers
if ((input%sigma_1_epsilon /= Zero) .and. (input%sigma_2_epsilon /= Zero)) then
  call WarningMessage(2,'SIG1 and SIG2 keywords are mutually exclusive')
  call Quit_OnUserError()
end if

if (input%sigma_1_epsilon > Zero) then
  sigma_p_epsilon = Input%sigma_1_epsilon
  sigma_p_exponent = 1
end if

if (input%sigma_2_epsilon > Zero) then
  sigma_p_epsilon = Input%sigma_2_epsilon
  sigma_p_exponent = 2
end if

do_real = real_shift > Zero
do_imag = imag_shift > Zero
do_sigp = sigma_p_epsilon > Zero
if ((do_real .and. (do_imag .or. do_sigp)) .or. (do_imag .and. do_sigp)) then
  call WarningMessage(2,'More than one intruder-state removal technique active: SHIFt/IMAGinary/SIG1/SIG2 are mutually exclusive!')
  call Quit_OnUserError()
end if

! RHS algorithm selection
if (Input%PRHS == '0') Input%PRHS = 'DEFAULT'
if (Input%PRHS == '1') Input%PRHS = 'OLD'
if (Input%PRHS == '2') Input%PRHS = 'NEW'
if (Input%PRHS == '3') Input%PRHS = 'DIRECT'  ! synonym of the DIREct keyword (undocumented)
if ((Input%PRHS /= 'DEFAULT') .and. (Input%PRHS /= 'OLD') .and. (Input%PRHS /= 'NEW') .and. (Input%PRHS /= 'DIRECT')) then
  call WarningMessage(1,'The selected PRHS is not supported. Going to use the default strategy.')
  Input%PRHS = 'DEFAULT'
end if
iParRHS = 1
#ifdef _MOLCAS_MPP_
! The RHS on-demand algorithm doesn't handle serial calculations
! because it's not adapted for use with regular Work arrays, only
! global arrays, and needs to be switched off (using rhsall instead)
RHSDIRECT = Is_Real_Par() .and. (Input%RHSD .or. (Input%PRHS == 'DIRECT'))
if (Is_Real_Par() .and. ((Input%PRHS == 'DEFAULT') .or. (Input%PRHS == 'NEW')) .and. (.not. RHSDIRECT)) iParRHS = 2
! maximum number of real values handled by a single GADGOP call
! compilers complain about huge(xxx)/RtoB with -Werror=integer-division
MAXBUF = (huge(1_MPIInt)-mod(int(huge(1_MPIInt),kind=iwp),RtoB))/RtoB
#else
RHSDIRECT = .false.
MAXBUF = (huge(RtoB)-mod(huge(RtoB),RtoB))/RtoB
#endif

! Cholesky: set defaults if it was not called during input
if (.not. (Input%ChoI .or. Input%Chol)) call Cho_caspt2_rdInp(.true.,iDummy)

! Initialize
RFpert = Input%RFPert

OUTFMT = 'DEFAULT'
G1SECIN = .false.
PRORB = Input%PrOrb
PRSD = Input%PrSD
NCASES = 13

JMS = Input%JMS
NLYROOT = Input%OnlyRoot
NLYGROUP = 0

DoCumulant = Input%DoCumulant
DMRG = Input%DMRG
CompressMPS = Input%CompressMPS

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

! then the arrays that are created will be:
! nstate = 3
! mstate(nstate): 3 4 5
! nGroup = 3
! nGroupState(nGroup): 1 1 1

! On the other hand, if the input is
! XMULtistate = 4 1 2 3 4

! then the arrays that are created will be:
! nstate = 4
! mstate(nstate): 1 2 3 4
! nGroup = 1
! nGroupState(nGroup): 4

NSTATE = 0
MSTATE(:) = 0
NGROUP = 0
NGROUPSTATE = 0
! This is the case for MS-CASPT2 and DW-CASPT2
if (Input%MULT) then
  if (Input%XMUL .or. Input%RMUL) then
    call WarningMessage(2,'Keyword MULT cannot be used with neither XMUL nor RMUL.')
    call Quit_OnUserError()
  end if
  ! Either the states were specified manually or the keyword "all"
  ! was used, so first we check the keyword all
  if (Input%AllMult) then
    NSTATE = NROOTS
    MSTATE(:) = IROOT
    NGROUP = NSTATE
    NGROUPSTATE(1:NGROUP) = 1
  else
    ! Save the states that need to be computed
    NGROUP = Input%nMultState
    MSTATE(1:NGROUP) = Input%MultGroup%A(1:NGROUP)
    NSTATE = NSTATE+NGROUP
    NGROUPSTATE(1:NGROUP) = 1
  end if
end if
IOFF = NSTATE

! This is the case for XMS-CASPT2 and XDW-CASPT2
if (Input%XMUL) then
  if (Input%MULT .or. Input%RMUL) then
    call WarningMessage(2,'Keyword XMUL cannot be used with neither MULT nor RMUL.')
    call Quit_OnUserError()
  end if
  ! This is a XDW-CASPT2 calculation. It is actually more similar to
  ! a MS-CASPT2 one since we need to put one state per group and thus
  ! have as many groups as states. Nevertheless, it makes more sense
  ! from the user point of view to ask for it through XMUL and DWMS
  if (Input%DWMS) then
    if (Input%AllXMult) then
      NSTATE = NROOTS
      MSTATE(:) = IROOT
      NGROUP = NSTATE
      NGROUPSTATE(1:NGROUP) = 1
    else
      NGROUP = Input%nXMulState
      MSTATE(1:NGROUP) = Input%XMulGroup%A(1:NGROUP)
      NSTATE = NSTATE+NGROUP
      NGROUPSTATE(1:NGROUP) = 1
    end if
    ! This is a XMS-CASPT2: one group with all the states
  else
    if (Input%AllXMult) then
      NSTATE = NROOTS
      MSTATE(:) = IROOT
      NGROUP = 1
      NGROUPSTATE(1) = NSTATE
    else
      NGROUP = 1
      NGROUPSTATE(NGROUP) = Input%nXMulState
      MSTATE(1:NGROUPSTATE(NGROUP)) = Input%XMulGroup%A(1:NGROUPSTATE(NGROUP))
      NSTATE = NSTATE+NGROUPSTATE(NGROUP)
    end if
  end if
end if

! This is the case for RMS-CASPT2
if (Input%RMUL) then
  if (Input%MULT .or. Input%XMUL) then
    call WarningMessage(2,'Keyword RMUL cannot be used with neither MULT nor XMUL.')
    call Quit_OnUserError()
  end if
  if (Input%AllRMult) then
    NSTATE = NROOTS
    MSTATE(:) = IROOT
    NGROUP = NSTATE
    NGROUPSTATE(1:NGROUP) = 1
  else
    NGROUP = Input%nRMulState
    MSTATE(1:NGROUP) = Input%RMulGroup%A(1:NGROUP)
    NSTATE = NSTATE+NGROUP
    NGROUPSTATE(1:NGROUP) = 1
  end if
end if
! After parsing mult, xmult or rmult, check that no two equal states where
! given in the input
if (Input%MULT .or. Input%XMUL .or. Input%RMUL) then
  do I=1,NSTATE
    do J=I+1,NSTATE
      if (MSTATE(I) == MSTATE(J)) then
        call WarningMessage(2,'The same root cannot be used twice in MULT/XMUL/RMUL blocks.')
        call Quit_OnUserError()
      end if
    end do
  end do
end if
! The LROOt keyword specifies a single root to be used. It should not be
! used together with either MULT or XMUL keywords.
if (Input%LROO) then
  if (Input%MULT .or. Input%XMUL .or. Input%RMUL) then
    call WarningMessage(2,'Keyword LROO cannot be used together with the MULT or XMUL keywords.')
    call Quit_OnUserError()
  end if
  NSTATE = 1
  MSTATE(1) = Input%SingleRoot
  NGROUP = 1
  NGROUPSTATE(1) = 1
end if
! if still nothing was selected we should default to compute all the
! roots that were part of the rasscf orbital optimization.
if (NSTATE == 0) then
  NSTATE = NROOTS
  MSTATE(:) = IROOT
  NGROUP = NSTATE
  NGROUPSTATE(1:NGROUP) = 1
end if
! Find the group number for OnlyRoot
if (NLYROOT /= 0) then
  IOFF = 0
  do IGROUP=1,NGROUP
    do I=1,NGROUPSTATE(IGROUP)
      if (IOFF+I == NLYROOT) NLYGROUP = IGROUP
    end do
    IOFF = IOFF+NGROUPSTATE(IGROUP)
  end do
end if
! Finally, some sanity checks.
if ((NSTATE <= 0) .or. (NSTATE > MXROOT)) then
  call WarningMessage(2,'Number of states is <0 or too large.')
  write(u6,'(a,i8)') ' NSTATE = ',NSTATE
  write(u6,*) ' Check usage of keywords MULT/XMUL/RMUL.'
  call Quit_OnUserError()
end if
! setup root to state translation
ROOT2STATE(:) = 0
do I=1,NSTATE
  ROOT2STATE(MSTATE(I)) = I
end do

! Relax root selection

iRlxRoot = -1
call Qpg_iScalar('NumGradRoot',Is_iRlxRoot_Set)
if (Is_iRlxRoot_Set) call Get_iScalar('NumGradRoot',iRlxRoot)
if (Input%RlxRoot > 0) iRlxRoot = Input%RlxRoot
if (iRlxRoot == -1) iRlxRoot = NSTATE
if (iRlxRoot > NSTATE) then
  if (IPRGLB >= TERSE) then
    call WarningMessage(1,'Too large iRlxRoot.')
    write(u6,*) ' Reset to NSTATE=',NSTATE
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
do iSym=1,nSym
  if (nCore(iSym) > nFro(iSym)) then
    nDiff = nCore(iSym)-nFro(iSym)
    nDiff = min(nDiff,nISh(iSym))
    nFro(iSym) = nFro(iSym)+nDiff
    nISh(iSym) = nISh(iSym)-nDiff
  end if
end do
! if a user specified the number of frozen orbitals explicitely, then
! that number overwrites the automatically chosen numbers, but it we
! warn the user and report what the computed number was.
if (Input%FROZ) then
  if (IPRGLB >= TERSE) call WarningMessage(1,'User changed nr of frozen orbitals.')
  do I=1,NSYM
    NFI = NFRO(I)+NISH(I)
    if (NFI < Input%nFro(I)) then
      call WarningMessage(2,'Too many frozen orbitals!')
      call Quit_OnUserError()
    else
      nFro(I) = Input%nFro(I)
    end if
    NISH(I) = NFI-nFro(I)
  end do
end if
! Set user-specified number of deleted orbitals.
if (Input%DELE) then
  do I=1,NSYM
    NSD = NSSH(I)+NDEL(I)
    if (NSD < Input%nDel(I)) then
      call WarningMessage(2,'Too many deleted orbitals!')
      call Quit_OnUserError()
    else
      NDEL(I) = Input%nDel(I)
    end if
    NSSH(I) = NSD-nDel(I)
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
  write(u6,*)
  write(u6,*) '  --------------------------------------'
  write(u6,*) '   Orbital-Free Embedding Calculation   '
  write(u6,*) '  --------------------------------------'
  write(u6,*)
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
IFMIX = .not. Input%NoMix
IFXMS = Input%XMUL
IFRMS = Input%RMUL
IFMSCOUP = (Input%MULT .or. IFXMS .or. IFRMS) .and. (.not. Input%NoMult)
if ((nState == 1) .and. (NLYROOT /= mState(1))) IFMSCOUP = .false.
IFDW = Input%DWMS
! Set type and exponent for DWMS
DWType = Input%DWType
if (IFDW) then
  if ((DWType <= 0) .or. (DWType > 3)) then
    call WarningMessage(2,'DWTYpe should be either 1, 2 or 3.')
    call Quit_OnUserError()
  end if
  zeta = Input%zeta
end if

! Choice? of preprocessing route
ORBIN = 'TRANSFOR'

! Expectation values of H0 needed or not? Spectral resolution?
BMATRIX = 'YES     '
BTRANS = 'YES     '
BSPECT = 'YES     '
! Overlap matrix needed? Linear dependence removal needed?
SMATRIX = 'YES     '
SDECOM = 'YES     '
! Presently, we always need the ON transformation/Lindep removal.
! Consistency of these demands:
if (BSPECT /= 'NO      ') BTRANS = 'YES     '
if (BTRANS /= 'NO      ') BMATRIX = 'YES     '
if (BTRANS /= 'NO      ') SDECOM = 'YES     '
if (SDECOM /= 'NO      ') SMATRIX = 'YES     '

#ifdef _DMRG_
if (DMRG) then
  ! just exit if somebody tries to do QD-CASPT2 with DMRG
  if (IFMSCOUP) then
    call WarningMessage(2,'Couplings with DMRG-CASPT2 not supported')
    call Quit_OnUserError()
  end if

  ! set OpenMolcas environment variables
  call getenv('Project',qcmaquis_param%project_name)
  call getenv('WorkDir',qcmaquis_param%workdir)

  ! save checkpoint paths for all states
  call mma_allocate(dmrg_file%qcmaquis_checkpoint_file,nstate)
  do i=1,nstate
    dmrg_file%qcmaquis_checkpoint_file(i) = trim(qcmaquis_param%workdir)//'/'//qcm_group_names(1)%states(mstate(i))
  end do

  ! initialize the interface using a checkpoint file
  if (iPrGlb >= VERBOSE) write(u6,*) 'PROCINP initializing QCMaquis DMRG interface...'
  call qcmaquis_interface_init_checkpoint(dmrg_file%qcmaquis_checkpoint_file(1))

  ! remove all measurements just to be sure
  call qcmaquis_interface_remove_param('MEASURE[1rdm]')
  call qcmaquis_interface_remove_param('MEASURE[2rdm]')
  call qcmaquis_interface_remove_param('MEASURE[3rdm]')
  call qcmaquis_interface_remove_param('MEASURE[4rdm]')
  call qcmaquis_interface_remove_param('MEASURE[1spdm]')
  call qcmaquis_interface_remove_param('MEASURE[ChemEntropy]')

  ! only set 1-rdm here and the other rdms later on in grpini
  call qcmaquis_interface_set_param('MEASURE[1rdm]','1')
  !call qcmaquis_interface_set_param('MEASURE[2rdm]','1')
  !call qcmaquis_interface_set_param('MEASURE[3rdm]','1')
  !call qcmaquis_interface_set_param('MEASURE[4rdm]','1')
end if
#endif
!***********************************************************************
!
! Thresholds
!
!***********************************************************************
if (Input%THRE) then
  if (IPRGLB >= TERSE) call WarningMessage(1,'User modified linear dependency thresholds!')
end if
THRSHN = Input%THRSHN
THRSHS = Input%THRSHS
THRCONV = Input%THRCONV
CITHR = Input%PrWF
THRENE = 50.0_wp
THROCC = 5.0e-4_wp
MAXIT = Input%MaxIter
DNMTHR = Input%DNMTHR
CMPTHR = Input%CMPTHR
CNTTHR = Input%CNTTHR

PT2Method = 'CASPT2'
CPT2Method = 'CASPT2'

Do_FIC = .true.
Do_SC = .false.
SC_prop = .false.
SC_amplitude = .false.
if (HZERO == 'DYALL') then
  MAXIT = 0
  Do_FIC = input%DOPC
  Do_SC = input%DOSC
  SC_prop = input%SCPROP .or. (.not. Do_FIC)
  SC_amplitude = SC_prop
  SC_thres = abs(input%SC_thres)
  if (NRAS1T+NRAS3T > 0) then
    call warningMessage(2,'NEVPT2 calculations with a RAS reference wavefunction are not supported')
    call quit_onUserError()
  end if
  if ((.not. Do_FIC) .and. (.not. Do_SC)) then
    call warningMessage(2,'Both NOPC and NOSC are specified? Cannot do the job...')
    call quit_onUserError()
  end if
  if (ipea_shift /= Zero) then
    ipea_shift = Zero
    if (IPRGLB >= TERSE) call warningMessage(1,'IPEA shift is automatically turned off for NEVPT2')
  end if
  if (do_real .or. do_imag .or. do_sigp) then
    if (IPRGLB >= TERSE) call warningMessage(1,'Real or imaginary shift or sigma regularization is being used only for PC-NEVPT2')
  end if
  ! The option XMULT forces to use the state-averaged Fock
  if (IFXMS) input%SADREF = .true.
  ! However, <Psi^0|Hact|Psi^0> is always diagonal, so skip the diagonalization process of the model states
  IFXMS = .false.
  IFRMS = .false.

  PT2Method = 'NEVPT2'
  CPT2Method = 'PC-NEVPT2'
  if (SC_prop) CPT2Method = 'SC-NEVPT2'
end if

!***********************************************************************
!
! Gradients
!
!***********************************************************************

! at the moment the calculation of analytic gradients has many
! technical restrictions and we need to make sure that we do not
! run into a unsupported combination of keywords, these are:
! 1. no symmetry
! 2. all CASSCF roots included in QD-CASPT2 gradients
! 3. QD-CASPT2 gradients only with CD/DF

call put_iScalar('mp2prpt',0)

! check if numerical gradients were requested in GATEWAY
call qpg_iScalar('DNG',DNG_available)
if (DNG_available) then
  call get_iScalar('DNG',iDNG)
  DNG = iDNG == 1
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
! I think GRDT keyword should be ignored for numerical gradient and last energy
do_lindep = .false.
if ((input%GRDT .or. input%NAC) .and. (SuperName(1:18) /= 'numerical_gradient') .and. (SuperName(1:11) /= 'last_energy')) then
  do_grad = Input%GRDT .or. input%NAC

  ! quit if both analytical and numerical gradients were explicitly requested
  if (DNG) then
    call warningMessage(2,'It seems that numerical gradients were requested in GATEWAY and analytical gradients in CASPT2. '// &
                        'Please choose only one of the two!')
    call quit_onUserError()
  end if

  ! only allow analytic gradients without symmetry
  if (nSym /= 1) then
    call warningMessage(2,'Analytic gradients only available without symmetry.')
    call quit_onUserError()
  end if

  if (ipea_shift /= Zero) do_lindep = .true.

  ! only allow analytic gradients either with nstate = nroots or with sadref
  if ((nState /= nRoots) .and. (.not. ifsadref)) then
    call warningMessage(2,'Analytic gradients available only if all CASSCF roots are included in the CASPT2 calculation or '// &
                        'with the SADRef keyword.')
    call quit_onUserError()
  end if

  ! QD-CASPT2 analytic gradients available only with DF or CD
  if (ifMSCoup .and. (.not. ifChol)) then
    call warningMessage(2,'MS-type analytic gradients available only with density fitting or Cholesky decomposition.')
    call quit_onUserError()
  end if

  ! CASPT2 analytic gradients with state-dependent density available only with DF or CD
  if ((.not. ifChol) .and. (.not. input%SADREF) .and. (nRoots /= 1)) then
    call warningMessage(2,'Analytic gradients with state-dependent density available only with density fitting or Cholesky '// &
                        'decomposition.')
    call quit_onUserError()
  end if
# ifdef _MOLCAS_MPP_
  ! No parallel without RI/CD
  if ((.not. ifChol) .and. (nProcs > 1)) then
    call warningMessage(2,'Analytic gradients without density fitting or Cholesky decomposition not available in parallel '// &
                        'executions.')
    call quit_onUserError()
  end if
# ifndef _GA_
  ! for the time being no gradients without GA
  ! Parallel CASPT2 gradient is implemented with some GA-specific subroutines
  ! partially because I can use OpenMolcas only for which GA is required.
  ! As long as OpenMolcas concerns, this should be no problem (at all)
  if (nProcs > 1) then
    call warningMessage(2,'Analytic gradients not available without GA installed. Install GA and link.')
    call quit_onUserError()
  end if
# endif
# endif

end if

! inside LAST_ENERGY we do not need analytic gradients
if (SuperName(1:11) == 'last_energy') DNG = .true.

! check if the calculation is inside a loop and make analytical
! gradients default in this case, unless the user specifically
! requested numerical gradients in GATEWAY
if (isStructure() == 1) then
  ! check the hard constraints first
  if ((.not. DNG) .and. (nSym == 1)) then
    do_grad = .true.

    ! check weaker constraints, if not met, revert to numerical gradients
    if (ifMSCoup .and. (.not. ifChol)) do_grad = .false.
    if ((ipea_shift /= Zero) .and. (.not. ifChol)) do_grad = .false.
    if ((nState /= nRoots) .and. (.not. ifsadref)) do_grad = .false.
  end if
end if

! compute full unrelaxed density for gradients
if (do_grad) ifDens = .true.

if (do_grad) then
  call put_iScalar('mp2prpt',2)
  do_nac = input%NAC

  !! If states to be computed are requested by ALASKA
  !! (if "@" presents in MCLR Root), always compute for these states
  call Qpg_cArray('MCLR Root',Found,I)
  if (Found) then
    call Get_cArray('MCLR Root',mstate1,16)
    if (mstate1 /= '****************') then
      if (index(mstate1,'@') /= 0) then
        read(mstate1,'(1X,I7,1X,I7)') iRoot1,iRoot2
        if (iRoot1 /= 0) do_nac = .true.
        if (iRoot1 == 0) then
          iRoot1 = iRoot2
          iRlxRoot = iRoot1
          do_nac = .false.
        end if
      end if
    end if
  end if

  !! If nothing is specified by ALASKA, use the states in &CASPT2
  if ((iRoot1 == 0) .and. (iRoot2 == 0)) then
    if ((input%iNACRoot1 == 0) .and. (input%iNACRoot2 == 0)) then
      iRoot1 = iRlxRoot
      iRoot2 = iRlxRoot
    else
      iRoot1 = input%iNACRoot1
      iRoot2 = input%iNACRoot2
    end if
  end if
end if

if (do_nac) then
  do_csf = input%CSF
  if (isStructure() == 1) do_csf = .false. !! omit the CSF term during (any) geometry optimizations
end if

IFSADREF = input%SADREF
IFDORTHO = input%DORTHO
if_invar = input%INVAR
if (ipea_shift /= Zero) if_invar = .false.
if_invaria = input%IAINVAR
if (sigma_p_epsilon /= Zero) if_invaria = .false. !! I'm not sure this is necessary, but it is needed for now
if (SC_prop) if_invaria = .false.
ConvInvar = input%ThrConvInvar

if ((ipea_shift /= Zero) .and. do_grad .and. (.not. IFDORTHO)) then
  call warningMessage(2,'Analytic gradients with IPEA shift must use the CORT or DORT option.')
  call quit_onUserError()
end if

#ifdef _MOLCAS_MPP_
if (do_grad .and. (sigma_p_epsilon /= Zero) .and. (nProcs > 1)) then
  call warningMessage(2,'Analytic gradients without the sigma^P regularization not available in parallel executions.')
  call quit_onUserError()
end if
#endif

if (do_grad .and. RF_On() .and. (.not. if_invar)) then
  if (IPRGLB >= TERSE) call warningMessage(1,'Analytic gradients with IPEA shift and PCM is not fully analytic.')
end if

!! Whether the Fock matrix (eigenvalues) is constructed with
!! the state-averaged density matrix or not.
!! The name of the variable is like state-specific DM,
!! but not necessarily state-specific. It is a matter of the
!! structure of DWGT(:,:) array or matrix.
!! DWGT is a matrix form for SS- and MS-CASPT2 with
!! state-specific DM, XDW-CASPT2, and RMS-CASPT2, while it is an
!! array for SS- and MS-CASPT2 with state-averaged DM (with SADREF
!! option) and XMS-CASPT2.
if (IFSADREF .or. (nRoots == 1) .or. (IFXMS .and. (.not. IFDW))) then
  if_SSDM = .false.
else
  if_SSDM = .true.
end if

!! Check if unequal-weighted MCSCF or not. Used only for gradients.
if (do_grad) then
  if (if_SSDM) if_equalW = .false.
  do I=2,nRoots
    if (Weight(1) /= Weight(I)) if_equalW = .false.
  end do
  if (.not. if_equalW) if_SSDM = .true.
end if

!! issue #448
if ((IFDENS .and. (.not. do_grad)) .and. (NRAS1T+NRAS3T > 0)) then
  call warningMessage(2,'DENS keyword cannot be combined with RAS.')
  call quit_onUserError()
end if

if (do_grad .and. (HZERO == 'DYALL')) then
  if ((.not. Do_FIC) .and. (.not. SC_prop)) then
    call warningMessage(2,'PC-NEVPT2 properties cannot be computed without PC-NEVPT2 energy. Remove the NOPC keyword.')
    call quit_onUserError()
  end if
  if ((.not. Do_SC) .and. SC_prop) then
    call warningMessage(2,'SC-NEVPT2 properties cannot be computed without SC-NEVPT2 energy. Remove the NOSC keyword.')
    call quit_onUserError()
  end if
end if

end subroutine procinp_caspt2
