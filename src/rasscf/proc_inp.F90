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
! Copyright (C) 2018, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine Proc_Inp(DSCF,lOPTO,iRc)

use Index_Functions, only: nTri_Elem
use fortran_strings, only: to_upper, operator(.in.)
use csfbas, only: CONF
use lucia_data, only: CFTP
use Fock_util_global, only: DoCholesky
use Cholesky, only: ChFracMem
use write_orbital_files, only: OrbFiles, write_orb_per_iter
use fcidump, only: DumpOnly
use fcidump_reorder, only: ReOrFlag, ReOrInp
use fciqmc, only: DoEmbdNECI, DoNECI, tGUGA_in, tNonDiagStochPT2, tPrepStochCASPT2
use fciqmc_read_RDM, only: MCM7, WRMA
use fciqmc_make_inp, only: definedet, diagshift, memoryfacspawn, nmCyc, pops_trial, RDMsampling, realspawncutoff, semi_stochastic, &
                           t_RDMsampling, Time, totalwalkers, trial_wavefunction
use CC_CI_mod, only: Do_CC_CI
use spin_correlation, only: orb_range_p, orb_range_q, same_orbs, tRootGrad
use orthonormalization, only: ON_scheme, ON_scheme_values
use casvb_global, only: ifvb
use KSDFT_Info, only: CoefR, CoefX
use OFembed, only: dFMD, Do_OFemb, KEonly, OFE_KSDFT, ThrFThaw, Xsigma
use CMS, only: CMSGiveOpt, CMSGuessFile, iCMSOpt
use UnixInfo, only: SuperName
use Lucia_Interface, only: Lucia_Util
use gugx, only: CIS, EXS, SGS
use gas_data, only: iDoGAS, IGSOCCX, NGAS, NGSSH
use Symmetry_info, only: Mul
use SplitCas_Data, only: DoSPlitCas, EnerSplit, fOrdSplit, GapSpli, iDimBlockA, lRootSplit, MxIterSplit, NumSplit, PerCSpli, &
                         PerSplit, ThrSplit
use PrintLevel, only: DEBUG, TERSE, VERBOSE
use output_ras, only: IPRGLB, IPRLOC
use general_data, only: CleanMask, CRPROJ, CRVec, INVEC, ISPIN, JOBIPH, JOBOLD, LOWDIN_ON, LUSTARTORB, MALTER, MAXALTER, NACTEL, &
                        NALTER, NASH, NBAS, NCONF, NCRVEC, NDEL, NDELT, NELEC3, NFRO, NFROT, NHOLE1, NISH, NORB, NRS1, NRS1T, &
                        NRS2, NRS2T, NRS3, NRS3T, NSEL, NSSH, NSYM, NTOT, NTOT1, NTOT2, NTOTSP, STARTORBFILE, STSYM, SXDAMP
use spinfo, only: DOBKAP, I2ELIMINATED_IN_GAS_MOLCAS, I_ELIMINATE_GAS_MOLCAS, IELIMINATED_IN_GAS_MOLCAS, IEXPAND_MOLCAS, &
                  IGSOCCX_MOLCAS, INOCALC_MOLCAS, IOCCPSPC, IPRCI_MOLCAS, IPT2_MOLCAS, ISAVE_EXP_MOLCAS, ISPEED, ISPIN_MOLCAS, &
                  ITMAX_MOLCAS, LSYM_MOLCAS, MS2, MS2_MOLCAS, N_2ELIMINATED_GAS_MOLCAS, N_ELIMINATED_GAS_MOLCAS, NACTEL_MOLCAS, &
                  NCSASM, NDET, NDTASM, NGAS_MOLCAS, NGASBK, NGSSH_MOLCAS, NROOTS_MOLCAS, NSYM_MOLCAS, POTNUC_MOLCAS, THRE_MOLCAS
use DWSol, only: DWSol_DWRO
use Molcas, only: LenIn, MxAct, MxGAS, MxOrb, MxRoot, MxSym
use RASDim, only: MxRef, MxTit
use input_ras, only: Key, LUInput, SetKey
use rasscf_global, only: CCI, CMSStartMat, CMSThreshold, CoreShift, DFTFOCK, DoBLOCKDMRG, DoFaro, DoFCIDump, ExFac, HFOCC, HFOcc, &
                         hRoots, iAlphaBeta, ICI, ICICH, iCIonly, iCIRFROOT, iCIRST, iCMSITERMAX, iCMSITERMin, iCMSP, iExpand, &
                         IfCRPR, iFORDE, InOCalc, iOrbOnly, iOrbTyp, iOrdEM, iOverWr, iPCMRoot, iPhName, iPR, iPT2, IRLXROOT, &
                         IROOT, iRotPsi, iSave_Exp, iSCF, iSPDEN, iSupSM, ITCORE, ITMAX, IXMSP, ixSym, iZRot, JCJ, kivo, KSDFT, &
                         kTight, l_CASDFT, LowMS, LROOTS, LvShft, MaxIt, MaxJt, MaxOrbOut, n_keep, NAC, NACPAR, NACPR2, NFR, NIN, &
                         NO2M, NonEq, NORBT, NQUNE, NROOTS, NSEC, NTOT3, NTOT4, OutFmt1, OutFmt2, PotNuc, PreThr, PreThr, ProThr, &
                         PrwThr, Purify, RFPert, S, SXSEL, ThFact, ThrE, ThrEn, ThrSX, ThrTE, Title, Weight
#ifdef _ENABLE_DICE_SHCI_
use rasscf_global, only: dice_eps1, dice_eps2, dice_iter, dice_restart, dice_sampleN, dice_stoc, diceOcc, nRef_dice
#endif
#ifdef _ENABLE_CHEMPS2_DMRG_
use rasscf_global, only: ChemPS2_BLB, ChemPS2_lRestart, ChemPS2_Noise, ChemPS2_Restart, Davidson_Tol, Do3RDM, Max_Canonical, &
                         Max_Sweep, MxDMRG
#endif
#ifdef _DMRG_
use qcmaquis_interface_cfg, only: dmrg_input, qcmaquis_param
use qcmaquis_interface, only: qcmaquis_interface_init, qcmaquis_interface_set_param, qcmaquis_interface_stdout, remove_comment
use active_space_solver_cfg, only: as_solver_inp_proc
use rasscf_global, only: DoDMRG, DoMCPDFTDMRG, DoNEVPT2Prep, MPSCompressM, Twordm_qcm
#ifdef _MOLCAS_MPP_
use Para_Info, only: mpp_nprocs, mpp_procid
#endif
#endif
#ifdef _HDF5_
use mh5, only: mh5_close_file, mh5_exists_attr, mh5_exists_dset, mh5_fetch_attr, mh5_fetch_dset, mh5_is_hdf5, mh5_open_file_r
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half, auToeV
use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp) :: DSCF, lOPTO
integer(kind=iwp) :: iRC
integer(kind=iwp) :: i, i_All, i1, i2, iad19, IADR19(15), iAlter, iChng1, iChng2, iDisk, iEnd, iErr, iGAS, iGrp, ii, ij, iJOB, &
                     iMAlter(MaxAlter,2), inporb_version, iod_save, iOffSet, iOrb, iOrbData, IPRGLB_IN, iPrLev, IPRLOC_IN(7), iR, &
                     iRC1, iRef, iReturn, is_in_group, IScratch(10), iStart, istatus, iSum, iSym, itu, j, jpcmroot, k, korb, kref, &
                     length, mBas, mCof, mConf, mm, mOrb, N, NA, NAO, NASHT, NBAS_L(8), NCHRG, nClean, nCof, NCRPROJ, NDEL_L(8), &
                     nDiff, NFRO_L(8), nGrp, NGSSH_HI, NGSSH_LO, NISH_L(8), NISHT, NISHT_old, nItems, nNUc, NORB_L(8), nOrbRoot, &
                     nOrbs, NRS1_L(8), NRS2_L(8), NRS3_L(8), NSSH_L(8), nSym_l, nT, nU, nW, start, step
real(kind=wp) :: dSum, dum1, dum2, dum3, Dummy(1), Eterna_1, Eterna_2, POTNUCDUMMY, PRO, SUHF, TEffNChrg, TotChrg
logical(kind=iwp) :: DBG, Exists, Langevin_On, lExists, PCM_On, RF_On, RlxRCheck, RunFile_Exists, SkipGUGA
character(len=(LenIn+8)*mxOrb) :: lJobH1
character(len=256) :: myTitle
character(len=180) :: Line
character(len=2*72) :: lJobH2
character(len=72) :: JobTit(mxTit), ReadStatus
character(len=50) :: ON_scheme_inp, uppercased
character(len=8) :: InfoLbl, MaxLab, NewJobIphName
integer(kind=iwp), allocatable :: iType(:), Stab(:), Temp1(:), Temp2(:), Temp3(:), UG2SG_X(:)
real(kind=wp), allocatable :: ENC(:), RF(:)
character(len=:), allocatable :: buffer
#ifdef _ENABLE_DICE_SHCI_
integer(kind=iwp) :: iref_dice
#endif
#ifdef _DMRG_
integer(kind=iwp) :: LRras2_dmrg(8), nr_lines
character(len=256) :: WorkDir
character(len=72) :: ProjectName
character(len=20) :: guess_dmrg
integer(kind=iwp), allocatable :: initial_occ(:,:)
#endif
#ifdef _HDF5_
integer(kind=iwp) :: lRoots_l, mh5id
character, allocatable :: typestring(:)
#endif
integer(kind=iwp), external :: Get_ExFac, IsFreeUnit, nToken
logical(kind=iwp), external :: Is_First_Iter
character(len=180), external :: Get_LN
#include "warnings.h"

!...Dongxia note for GAS:
!   No changing about read in orbital information from INPORB yet.

DoFaro = .false.

#ifdef _DMRG_
! Leon: Prepare 4-RDM calculations for (CD)-DMRG-NEVPT2 at the end of the calculation
DoNEVPT2Prep = .false.
! If this is set to 0, MPS compression is disabled
MPSCompressM = 0
#endif
! NN.14 Block DMRG flag
DoBlockDMRG = .false.
#ifdef _ENABLE_CHEMPS2_DMRG_
! Quan.16: CheMPS2 default flags
chemps2_restart = .false.
chemps2_lrestart = 0
davidson_tol = 1.0e-7_wp
chemps2_blb = 0.5e-2_wp
max_sweep = 8
chemps2_noise = 0.05_wp
max_canonical = max_sweep*5
#endif
! Init HFOCC array containing user defined occupancies for the active orbitals.
! This array is used by DMRG codes (Block as well as CheMPS2).
! Terefore I took it out of any ifdef preprocessing flag.

hfocc(:) = 0

#ifdef _ENABLE_DICE_SHCI_
dice_stoc = .false.
nref_dice = 1
dice_eps1 = 1.0e-4_wp
dice_eps2 = 1.0e-5_wp
dice_sampleN = 200
dice_iter = 20
dice_restart = .false.
#endif

! SplitCAS related variables declaration  (GLMJ)
DoSplitCAS = .false.
NumSplit = .false.
EnerSplit = .false.
PerSplit = .false.
FOrdSplit = .false.
! BK type of approximation (GLMJ)
DoBKAP = .false.

! ======================================================================
!   QCMaquis flags
! ======================================================================
dofcidump = .false.

! ======================================================================

! GAS flag, means the INPUT was GAS
iDoGas = .false.

! The compiler thinks NASHT could be undefined later (after 100)
NASHT = 0

DBG = .false.
NAlter = 0
iRc = _RC_ALL_IS_WELL_

! Note: During process of keywords, there is either a call to ChkIfKey
! once any input data for that keyword has been processed, to check that
! there are no unrecognized keywords (misspelled? syntactic error?)
! following; or if the keyword should not be followed by any data,
! then for the same purpose, a call of the kind
!call SetPos(LUInput,'ATOM',Line,iRc)
!call ChkIfKey()
! This is really the only safe way to make warnings for such events
! the way the input is parsed and processed for the moment.
! Start by checking for extraneous input following any of these
! keywords (misspelled or nonsyntactic usage)
if (Key('CORE')) then
  call SetPos(LUInput,'CORE',Line,iRc)
  call ChkIfKey()
end if
if (Key('INPO')) then
  call SetPos(LUInput,'INPO',Line,iRc)
  call ChkIfKey()
end if
if (Key('JOBI')) then
  call SetPos(LUInput,'JOBI',Line,iRc)
  call ChkIfKey()
end if
if (Key('LUMO')) then
  call SetPos(LUInput,'LUMO',Line,iRc)
  call ChkIfKey()
end if
if (Key('CIRE')) then
  call SetPos(LUInput,'CIRE',Line,iRc)
  call ChkIfKey()
end if
! They had to be tested here, since the key flags may be modified later.

! ======================================================================
! The outdated INPORB keyword now means just the same as LUMORB
! Will be deleted soon.
call SetKey('LUMO',Key('LUMO') .or. Key('INPO'))
! ======================================================================

! How was the program called?
!PAM 2009  For some particular types of calculations, the input is to be
! disregarded or overridden, as follows (Roland, esp. numerical diff):
if (Key('EXPE')) then
  call SetPos(LUInput,'EXPE',Line,iRc)
  call ChkIfKey()
else
  IfVB = 0
  if (SuperName(1:6) == 'rasscf') then
    ! For geometry optimizations use the old CI coefficients.
    if (.not. Is_First_Iter()) then
      call SetKey('CIRE',.true.)
      call SetKey('FILE',.false.)
    end if
  else if (SuperName(1:5) == 'casvb') then
    IfVB = 2
  else if (SuperName(1:6) == 'loprop') then
    call SetKey('CIRE',.true.)
    call SetKey('FILE',.false.)
  else if (SuperName(1:11) == 'last_energy') then
    call SetKey('CIRE',.true.)
    call SetKey('FILE',.false.)
  else if (SuperName(1:18) == 'numerical_gradient') then
    call SetKey('CIRE',.true.)
    call SetKey('FILE',.false.)
  end if
end if
!PAM2009 Also, CIRESTART would normally also imply that orbitals are
!        to be taken from JOBIPH or JOBOLD:
if (.not. Key('EXPE')) then
  if (Key('CIRE')) then
    call SetKey('LUMO',.false.)
    call SetKey('INPO',.false.)
    call SetKey('JOBI',.true.)
  end if
end if

! Check PRINT command
if (Key('PRIN')) then
  call SetPos(LUInput,'PRIN',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) then
    call Error(1)
    return
  end if
  ReadStatus = ' Failure reading after PRINT keyword.'
  read(LUInput,*,iostat=istatus) j,IPRGLB_IN,(IPRLOC_IN(I),I=1,7)
  if (istatus < 0) then
    call Error(2)
    return
  else if (istatus > 0) then
    call Error(3)
    return
  end if
  ReadStatus = ' O.K. reading info of PRINT keyword.'
  call SetPrLev(IPRGLB_IN,IPRLOC_IN)
  if (IPRLOC(1) >= DEBUG) then
    write(u6,*) ' PROC_INP: Print levels have been set:'
    write(u6,*) '    Global print level IPRGLB=',IPRGLB
    write(u6,*) '    Local print levels by section:'
    write(u6,*) '             Input section, IPRLOC(1)=',IPRLOC(1)
    write(u6,*) '    Transformation section, IPRLOC(2)=',IPRLOC(2)
    write(u6,*) '                CI section, IPRLOC(3)=',IPRLOC(3)
    write(u6,*) '          Super-CI section, IPRLOC(4)=',IPRLOC(4)
    write(u6,*) '            Output section, IPRLOC(6)=',IPRLOC(6)
    write(u6,*) '          Property section, IPRLOC(7)=',IPRLOC(7)
  end if
  call ChkIfKey()
end if

! Local print level in this routine:
IPRLEV = IPRLOC(1)
! Short, for use with IF-statements at debugging print level:
DBG = DBG .or. (IPRLEV >= DEBUG)

if (DBG) write(u6,*) ' Trace of input processing in PROC_INP:'

! ======================================================================
! Check IfVB flag
! if (IFVB == 2) -- Then bypass a lot of processing:
if (DBG) write(u6,*) ' Valence-Bond flag IfVB=',IfVB
if (ifvb == 2) then
  if (JOBIPH <= 0) then
    JOBIPH = IsFreeUnit(15)
    call daname(JOBIPH,'JOBIPH')
  end if
  ! Special input routine:
  if (DBG) write(u6,*) ' IfVB=2, so Call Readin_VB'
  call Readin_vb()
  if (DBG) then
    write(u6,*) ' Back from Readin_VB'
    write(u6,*) ' Bypass usual input processing!'
  end if
else

  ! ====================================================================
  if (DBG) write(u6,*) ' Check if VB keyword was given.'
  !---  Process vb   command
  if (Key('VB')) then
    if (DBG) write(u6,*) ' Yes it was!'
    call SetPos(LUInput,'VB  ',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    if (DBG) write(u6,*) ' so Call CVBInp_CVB.'
    call cvbinp_rvb(1,LUInput)
    if (DoCholesky) then
      call WarningMessage(2,'CASVB cannot do Cholesky or RI/DF.')
      call Quit(_RC_INPUT_ERROR_)
    end if
    if (DBG) write(u6,*) ' Set IfVB=1.'
    ifvb = 1
    if (DBG) write(u6,*) ' Asked for CASVB calc.'
  end if

  ! ====================================================================
  !---  process TITLE    command
  if (Key('TITL')) then
    call SetPos(LUInput,'TITL',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    Title(1) = trim(Get_Ln(LUInput))
    if (DBG) then
      write(u6,*) ' PROC_INP: Title line:'
      write(u6,*) Title(1)
    end if
    call ChkIfKey()
  end if

  ! ====================================================================
  ! Check if there is any runfile:
  call F_Inquire('RUNFILE',RunFile_Exists)
  if (DBG) write(u6,*) ' Inquire about RUNFILE.'
  if (RunFile_Exists) then
    if (DBG) write(u6,*) ' Yes, there is one.'
    NSYM = 0
    call qpg_iScalar('nSym',lExists)
    if (lExists) then
      call Get_iScalar('nSym',nSym)
      call Get_iArray('nBas',nBas,nSym)
      if (DBG) then
        write(u6,*) ' The following information exists on runfile:'
        write(u6,*) ' Nr of symmetries, NSYM:',NSYM
        write(u6,*) ' Nr of basis functions/symmetry:'
        write(u6,'(1x,8I5)') (NBAS(I),I=1,NSYM)
      end if
    else
      call WarningMessage(2,'No symmetry info on runfile.')
      write(u6,*) ' There seems to be no information about symmetry'
      write(u6,*) ' on the runfile! This is an unexpected error.'
      call Quit(_RC_IO_ERROR_READ_)
    end if
  else
    call WarningMessage(2,'Cannot find runfile.')
    write(u6,*) ' PROC_INP: Cannot find RUNFILE. This is an unexpected error.'
    call Quit(_RC_IO_ERROR_READ_)
  end if

  ! PAM03: Local variable IORBDATA=0 as long as no explicit orbital specifications have
  ! been set by the user input.
  ! iOrbData=0: no orbital space data is specified
  !          1: user input has defined orbital spaces
  !          2: specifications from runfile, or JOBOLD or JOBIPH file
  !          3: specifications from orbital file
  iOrbData = 0
  INVEC = 0
  ! INVEC=0, no source for orbitals (yet)
  !       1, CORE command: compute orbitals from scratch.
  !       2, read from starting orbitals file in INPORB format.
  !       3, take from JOBOLD, or JOBIPH file
  ! RASSCF will read CMO arrays by READVC call after return from PROC_INP
  ! Inside readvc, if INVEC is still 0, it will be changed to 4, 5, or 6:
  !       4, take from an HDF5 file
  !       5, take rasscf orbitals from runfile
  !       6, take scf orbitals from runfile
  !       7, take guessorb orbitals from runfile
  ! iOverWr: If orbital coefficient arrays (CMO:s) are read from INPORB
  ! file, then iOverWr=0 means that the arrays are reordered by typeindex.
  ! iOverWr=1 means do not reorder by typeindex.
  iOverWr = 0

  ! Where should orbital be read (if at all)?
  if (DBG) write(u6,*) ' Where to read orbitals? '
  StartOrbFile = 'INPORB'
  if (Key('FILE')) then
    if (DBG) write(u6,*) ' Reading file name for start orbitals.'
    call SetPos(LUInput,'FILE',Line,iRc)
    Line = Get_Ln(LUInput)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    call ChkIfKey()
    if (DBG) then
      write(u6,*) ' Calling fileorb with filename='
      write(u6,*) Line
    end if
    call fileorb(Line,StartOrbFile)
#   ifdef _HDF5_
    if (mh5_is_hdf5(StartOrbFile)) then
      call SetKey('LUMO',.false.)
      call SetKey('TYPE',.false.)
      call SetKey('H5OR',.true.)
      if (Key('CIRE')) then
        write(u6,*) ' CIRE keyword given, but FILE is an HDF5 format.'
        write(u6,*) ' I will switch off CIRE and use H5CI instead!'
        call SetKey('H5CI',.true.)
        call SetKey('CIRE',.false.)
        call SetKey('JOBI',.false.)
      end if
    else
      if (Key('CIRE')) then
        write(u6,*) ' CIRE keyword given, FILE will be ignored!'
      else
        call SetKey('LUMO',.true.)
      end if
    end if
#   else
    if (Key('CIRE')) then
      write(u6,*) ' CIRE keyword given, FILE will be ignored!'
    else
      call SetKey('LUMO',.true.)
    end if
#   endif
  end if
  if (DBG) then
    write(u6,*) ' StartOrbFile='
    write(u6,*) StartOrbFile
  end if

  ! ====================================================================
  ! The JOBIPH file, following decisions from the Torre Normanna labour camp:
  ! Default, the file name is 'JOBIPH'.
  ! However, if keyword IPHNAME was used, then the name was given in input.
  ! Also, if instead the keyword NEWIPH was given, then a new name will be
  ! chosen as the first not-already-used name in the sequence
  ! 'JOBIPH', 'JOBIPH01', 'JOBIPH02', etc.
  if (DBG) write(u6,*) ' Present name of JOBIPH file is ',IPHNAME
  !---  Process NEWIPH command (P A Malmqvist Sep 06)
  if (Key('NEWI')) then
    if (DBG) write(u6,*) ' A fresh JOBIPH file will be used.'
    IPHNAME = 'ToBeFoun'
  end if
  !---  process IPHNAME command (P A Malmqvist Sep 06)
  if (Key('IPHN')) then
    if (DBG) write(u6,*) ' Reading file name for JOBIPH file.'
    call SetPos(LUInput,'IPHN',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading IPHNAME string.'
    read(LUInput,*,iostat=istatus) IPHNAME
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading IPHNAME string.'
    call UpCase(IPHNAME)
  end if
  ! The JOBIPH file:
  if (IPHNAME == 'ToBeFound') then
    ! Choose a new jobiph name.
    call F_Inquire('JOBIPH',lExists)
    if (.not. lExists) then
      IPHNAME = 'JOBIPH'
    else
      do I=1,99
        write(NewJobIphName,'(a6,i2.2)') 'JOBIPH',I
        call F_Inquire(NewJobIphName,lExists)
        ! Choose first non-existent file
        if (.not. lExists) exit
      end do
      if (I > 99) then
        write(u6,*)
        write(u6,*) '******************************************'
        write(u6,*) ' Sorry, all possible JOBIPH names in this '
        write(u6,*) ' directory already in use. Program stops. '
        write(u6,*) '******************************************'
        call Abend()
      end if
      IPHNAME = NewJobIphName
    end if
  end if
  if (DBG) then
    write(u6,*) ' Name of JOBIPH file is'
    write(u6,*) IPHNAME
  end if
  ! Finally, we know the name. Open jobiph file. If another file is already
  ! opened with this identifier close it.
  if (JOBIPH > 0) then
    call DaClos(JOBIPH)
    JOBIPH = -1
  end if
  JOBIPH = IsFreeUnit(15)
  call DANAME(JOBIPH,IPHNAME)

  ! ====================================================================
  ! If orbital files should be produced, several keywords are relevant:
  !---  Process OUTPRINT command
  if (Key('OUTP')) then
    if (DBG) write(u6,*) ' OUTPRINT command was given:'
    call SetPos(LUInput,'OUTP',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading OUTPRINT line.'
    read(LUInput,*,iostat=istatus) Line
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading OUTPRINT line.'
    call UpCase(Line)
    if (index(Line,'VERY') /= 0) OutFmt1 = 'NOTHING'
    if (index(Line,'BRIE') /= 0) OutFmt1 = 'FEW'
    if (index(Line,'LONG') /= 0) OutFmt1 = 'ALL'
    if (index(Line,'NOTH') /= 0) OutFmt1 = 'NOTHING'
    if (index(Line,'FEW') /= 0) OutFmt1 = 'FEW'
    if (index(Line,'NOCO') /= 0) OutFmt1 = 'NOCORE'
    if (index(Line,'ALL') /= 0) OutFmt1 = 'ALL'
    if (index(Line,'COMP') /= 0) OutFmt2 = 'COMPACT'
    if (index(Line,'FULL') /= 0) OutFmt2 = 'FULL'
    if ((OutFmt1 == 'DEFAULT') .and. (OutFmt2 == 'DEFAULT')) then
      call WarningMessage(1,'Error in ''OUTP'' command?')
      write(u6,*) ' Input line is:'
      write(u6,*) Line
      write(u6,*) ' Did not understand ''OUTP'' input. Ignored.'
    end if
  end if
  !---  Process PROR (Print levels for orbitals) command
  if (Key('PROR')) then
    if (DBG) write(u6,*) ' PRORB command was given:'
    call SetPos(LUInput,'PROR',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading PRORB input.'
    read(LUInput,*,iostat=istatus) PRETHR,PRO
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading PRORB input.'
    PROTHR = max(Zero,PRO)
  end if
  !---  Process ORBL keyword: Orbital listing
  if (Key('ORBL')) then
    if (DBG) write(u6,*) ' ORBL (Orbital listing):'
    call SetPos(LUInput,'ORBL',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading ORBL input.'
    read(LUInput,*,iostat=istatus) Line
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading ORBL input.'
    call UpCase(Line)
    if (index(Line,'VERY') /= 0) OutFmt1 = 'NOTHING '
    if (index(Line,'BRIE') /= 0) OutFmt1 = 'FEW     '
    if (index(Line,'LONG') /= 0) OutFmt1 = 'ALL     '
    if (index(Line,'NOTH') /= 0) OutFmt1 = 'NOTHING '
    if (index(Line,'FEW') /= 0) OutFmt1 = 'FEW     '
    if (index(Line,'NOCO') /= 0) OutFmt1 = 'NOCORE  '
    if (index(Line,'ALL') /= 0) OutFmt1 = 'ALL     '
    if (OutFmt1 == 'DEFAULT ') then
      call WarningMessage(1,'Error in ''ORBL'' command?')
      write(u6,*) ' Input line is:'
      write(u6,*) Line
      write(u6,*) ' Did not understand ''OUTL'' input. Ignored.'
    end if
  end if
  !---  Process ORBA keyword: Orbital Appearance
  if (Key('ORBA')) then
    if (DBG) write(u6,*) ' ORBA (Orbital appearance):'
    call SetPos(LUInput,'ORBA',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading ORBA input.'
    read(LUInput,*,iostat=istatus) Line
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading ORBA input.'
    call UpCase(Line)
    if (index(Line,'COMP') /= 0) OutFmt2 = 'COMPACT '
    if (index(Line,'FULL') /= 0) OutFmt2 = 'FULL    '
    if (OutFmt2 == 'DEFAULT ') then
      call WarningMessage(1,'Error in ''OUTA'' command?')
      write(u6,*) ' Input line is:'
      write(u6,*) Line
      write(u6,*) ' Did not understand ''OUTA'' input. Ignored.'
    end if
  end if
  !---  Process MAXO keyword: Max nr of state-specific orbital files produced
  if (Key('MAXO')) then
    if (DBG) write(u6,*) ' MAXORB command:'
    call SetPos(LUInput,'MAXO',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading MAXO input.'
    read(LUInput,*,iostat=istatus) MAXORBOUT
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading MAXO input.'
  end if
  if (DBG) then
    write(u6,*) ' Orbital print levels are'
    write(u6,*) '   for energy, PRETHR=',PRETHR
    write(u6,*) '   for occup , PROTHR=',PROTHR
    write(u6,*) ' Orbital listing flag is: ',OutFmt1
    write(u6,*) ' Orbital appearance flag: ',OutFmt2
    write(u6,*) ' Max nr of state-specific orbital files is ',MAXORBOUT
  end if
  ! ====================================================================
  !---  Process OUTORBITALS command (Which kind of orbitals)
  if (Key('OUTO')) then
    call SetPos(LUInput,'OUTO',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading OUTO input.'
    read(LUInput,*,iostat=istatus) Line
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading OUTO input.'
    call UpCase(Line)
    iOrbTyp = 0
    if (index(Line,'AVER') /= 0) iOrbTyp = 1
    if (index(Line,'CANO') /= 0) iOrbTyp = 2
    if (index(Line,'NATU') /= 0) iOrbTyp = 3
    if (index(Line,'SPIN') /= 0) iOrbTyp = 4
    if (iOrbTyp == 0) then
      write(u6,*) ' The line after keyword ''OUTORBITALS'' is'
      write(u6,*) ' not understood. That line begins:'
      write(u6,'(1x,a60)') line(1:60)
      write(u6,*) ' This input is IGNORED.'
    end if
    if (iOrbTyp == 2) IPT2 = 1
    if ((iOrbTyp == 3) .or. (iOrbTyp == 4)) then
      ReadStatus = ' Failure reading nOrbRoot after OUTO keyword.'
      read(LUInput,*,iostat=istatus) nOrbRoot
      if (istatus < 0) then
        call Error(2)
        return
      else if (istatus > 0) then
        call Error(3)
        return
      end if
      ReadStatus = ' O.K. after reading nOrbRoot after OUTO keyword.'
    end if
  end if
  if (DBG) then
    write(u6,*) ' OUTORBITALS command specified orbital type ',iOrbTyp
    if (iOrbTyp == 1) write(u6,*) ' Meaning: Average'
    if (iOrbTyp == 2) write(u6,*) ' Meaning: Canonical'
    if (iOrbTyp == 3) write(u6,*) ' Meaning: Natural'
    if (iOrbTyp == 4) write(u6,*) ' Meaning: Spin orbitals.'
    if ((iOrbTyp == 3) .or. (iOrbTyp == 4)) write(u6,*) ' Max state for printing this orbital type ',nOrbRoot
  end if
  !---  Process ORDER command (SVC Feb 06)
  if (Key('ORDE')) then
    if (DBG) write(u6,*) ' ORDER command was used.'
    call SetPos(LUInput,'ORDE',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure after reading ORDER keyword.'
    read(LUInput,*,iostat=istatus) IFORDE
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading ORDER keyword.'
    IORDEM = 1
    if (DBG) write(u6,*) ' IFORDE, IORDEM=',IFORDE,IORDEM
  end if
  !---  Process PRSP command
  if (Key('PRSP')) then
    if (DBG) write(u6,*) ' PRSPIN command was used.'
    ISPDEN = 1
  end if
  !---  Process IVO command
  if (Key('IVO')) then
    if (DBG) write(u6,*) ' IVO command was used.'
    kIvo = .true.
  end if
  ! ====================================================================
  !  If ORBONLY keyword was used, then the JOBIPH file should be used
  ! only to produce orbital files, then the program stops.
  !---  Process ORBO command -( new! generate orbitals from jobiph, only)
  if (Key('ORBO')) then
    iOrbOnly = 1
    call OrbFiles(JOBIPH,IPRLEV)
    if (JOBIPH > 0) then
      call DaClos(JOBIPH)
      JOBIPH = -1
    end if
    ! Nothing more to be done, so return.
    iReturn = _RC_ALL_IS_WELL_
    call xQuit(iReturn)
  end if

  !---  Process ALTEr command (G. Ghigo Sep 03)
  if (Key('ALTE')) then
    if (DBG) write(u6,*) ' ALTER command has been used.'
    call SetPos(LUInput,'ALTE',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure after reading ALTER keyword.'
    read(LUInput,*,iostat=istatus) NAlter
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading ALTER keyword.'
    if (NAlter > MaxAlter) then
      write(MaxLab,*) MaxAlter
      MaxLab = adjustl(MaxLab)
      write(u6,*)
      call WarningMessage(2,'Alter: too many orbital pairs.')
      write(u6,*) ' ************* ERROR **************'
      write(u6,*) ' ALTEr: Too many pairs of orbitals '
      write(u6,*) ' to exchange (max '//trim(MaxLab)//').'
      write(u6,*) ' **********************************'
      call Abend()
    end if
    do iAlter=1,NAlter
      ReadStatus = ' Failure reading data after ALTER keyword.'
      read(LUInput,*,iostat=istatus) (MAlter(iAlter,i),i=1,3)
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
      ReadStatus = ' O.K. after reading data after ALTER keyword.'
    end do
    ! (SVC) get absolute orbital values for the alterations so that
    ! iMAlter is symmetry independent
    if (DBG) write(u6,*) ' ''Absolute'' iMAlter indices:'
    do iAlter=1,NAlter
      iEnd = 0
      iStart = 1
      do iSym=1,MAlter(iAlter,1)
        iStart = iEnd+1
        iEnd = iEnd+nBas(iSym)
      end do
      iMAlter(iAlter,1) = MAlter(iAlter,2)+iStart-1
      iMAlter(iAlter,2) = MAlter(iAlter,3)+iStart-1
      if (DBG) write(u6,'(1x,2I5)') iMAlter(iAlter,1),iMAlter(iAlter,2)
    end do
  end if
  !---  Process ATOM command (P A Malmqvist Apr 07)
  if (Key('ATOM')) then
    PURIFY = 'ATOM    '
    ISUPSM = 1
    call SetPos(LUInput,'ATOM',Line,iRc)
    call ChkIfKey()
  end if
  !---  Process LINEAR command (P A Malmqvist Apr 05)
  if (Key('LINE')) then
    PURIFY = 'LINEAR'
    ISUPSM = 1
    call SetPos(LUInput,'LINE',Line,iRc)
    call ChkIfKey()
  end if
  if (DBG) write(u6,*) ' Purify=',PURIFY

  !---  process KSDF command
  if (DBG) write(u6,*) ' Check if FUNCtional was requested.'
  if (Key('FUNC')) then
    if (DBG) write(u6,*) ' FUNC command was given.'
    DFTFOCK = 'CAS '
    call SetPos(LUInput,'FUNC',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    read(LUInput,*,iostat=istatus) Line
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    call UpCase(Line)
    if (Line(1:4) == 'ROKS') DFTFOCK = 'ROKS'
    if (Line(1:6) == 'CASDFT') DFTFOCK = 'DIFF'
    read(LUInput,*,iostat=istatus) Line
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    KSDFT = Line(1:80)
    call UpCase(KSDFT)
    l_casdft = (KSDFT(1:2) == 'T:') .or. (KSDFT(1:3) == 'FT:')
    if (.not. l_casdft) then
      call Error(3)
      return
    end if
    if ((IPRLOC(1) >= DEBUG) .and. l_casdft) write(u6,*) ' MCPDFT with functional:',KSDFT
    ExFac = Get_ExFac(KSDFT)
    !---  Process DFCF command
    if (Key('DFCF')) then
      if (DBG) write(u6,*) ' DFCF (dens. func. coeff) command was used.'
      call SetPos(LUInput,'DFCF',Line,iRc)
      if (iRc /= _RC_ALL_IS_WELL_) then
        call Error(1)
        return
      end if
      ReadStatus = ' Failure reading DF coeff after DFCF keyword.'
      read(LUInput,*,iostat=istatus) CoefX,CoefR
      if (istatus < 0) then
        call Error(2)
        return
      else if (istatus > 0) then
        call Error(3)
        return
      end if
      ReadStatus = ' O.K. after reading DF coeff after DFCF keyword.'
      if (DBG) then
        write(u6,*) ' Density functional exchange coeff, CoefX=',CoefX
        write(u6,*) ' Density functional correlation coeff, CoefRE=',CoefR
      end if
      call ChkIfKey()
    end if

    call ChkIfKey()
  end if
  !---  Process CION command
  if (DBG) write(u6,*) ' Check if CIONLY case.'
  if (Key('CION')) then
    if (DBG) write(u6,*) ' CIONLY keyword was used.'
    iCIonly = 1
    call SetPos(LUInput,'CION',Line,iRc)
    call ChkIfKey()
  end if
  !---  Process ROST command
  if (DBG) write(u6,*) ' Check if ROtSTate case.'
  if (Key('ROST')) then
    if (DBG) write(u6,*) ' ROtSTate keyword was used.'
    iRotPsi = 1
    call SetPos(LUInput,'ROST',Line,iRc)
    call ChkIfKey()
  end if
  !---  Process XMSI command
  if (DBG) write(u6,*) ' Check if XMSI case.'
  if (Key('XMSI')) then
    if (DBG) write(u6,*) ' XMSI keyword was used.'
    iRotPsi = 1
    IXMSP = 1
    call SetPos(LUInput,'XMSI',Line,iRc)
    call ChkIfKey()
  end if
  !---  Process CMSI command
  if (DBG) write(u6,*) ' Check if CMSI case.'
  if (Key('CMSI')) then
    if (DBG) write(u6,*) ' CMSI keyword was used.'
    iRotPsi = 1
    ICMSP = 1
    call SetPos(LUInput,'CMSI',Line,iRc)
    call ChkIfKey()
  end if
  !---  Process CMSS command
  CMSStartMat = 'XMS'
  if (Key('CMSS') .and. (iCMSP == 1)) then
    if (DBG) write(u6,*) ' Reading CMS initial rotation matrix'
    call SetPos(LUInput,'CMSS',Line,iRc)
    Line = Get_Ln(LUInput)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    call ChkIfKey()
    if (DBG) then
      write(u6,*) ' Reading CMS starting rotation matrix from'
      write(u6,*) trim(Line)
    end if
    if (.not. (Line == 'XMS')) then
      CMSGuessFile = trim(Line)
      CMSStartMat = CMSGuessFile
      call F_Inquire(trim(CMSStartMat),lExists)
      if (.not. lExists) then
        write(u6,'(6X,A,A)') trim(CMSStartMat),' is not found. Use XMS intermediate states as initial guess.'
        CMSStartMat = 'XMS'
      end if
      !call fileorb(Line,CMSStartMat)
    end if
  end if
  !---  Process CMSO command
  if (Key('CMSO') .and. (iCMSP == 1)) then
    if (DBG) write(u6,*) 'Inputting CMS optimization option'
    call SetPos(LUInput,'CMSO',Line,iRc)
    Line = Get_Ln(LUInput)
    call Upcase(Line)
    if (Line(1:4) == 'NEWT') then
      iCMSOpt = 1
    else if (Line(1:4) == 'JACO') then
      iCMSOpt = 2
    else
      ReadStatus = 'Wrong value assigned to keyword CMSO'
      call Error(3)
      return
    end if
    CMSGiveOpt = .true.
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    if (DBG) write(u6,*) ' CMS Optimization Option',iCMSOpt
    call ChkIfKey()
  end if
  !---  Process CMMA command
  if (Key('CMMA')) then
    if (DBG) write(u6,*) ' CMS Max Cylces keyword was given.'
    call SetPos(LUInput,'CMMA',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading data following CMMA keyword.'
    read(LUInput,*,iostat=istatus) ICMSIterMax
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. reading data following CMMA keyword.'
    if (DBG) write(u6,*) ' Max nr of CMS cycles',ICMSIterMax
    call ChkIfKey()
  end if
  !---  Process CMMI command
  if (Key('CMMI')) then
    if (DBG) write(u6,*) ' CMS Min Cylces keyword was given.'
    call SetPos(LUInput,'CMMI',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading data following CMMI keyword.'
    read(LUInput,*,iostat=istatus) ICMSIterMin
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. reading data following CMMI keyword.'
    if (DBG) write(u6,*) ' Min nr of CMS cycles',ICMSIterMin
    call ChkIfKey()
  end if
  !---  Process CMTH command
  if (Key('CMTH')) then
    if (DBG) write(u6,*) ' CMS Threshold keyword was given.'
    call SetPos(LUInput,'CMTH',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading data following CMTH keyword.'
    read(LUInput,*,iostat=istatus) CMSThreshold
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. reading data following CMTH keyword.'
    if (DBG) write(u6,*) ' CMS threshold',CMSThreshold
    call ChkIfKey()
  end if
  !---  Process RFPE command
  if (Key('RFPE')) then
    if (DBG) write(u6,*) ' RFPERT (Response Field Perturbation)'
    RFpert = .true.
    call SetPos(LUInput,'RFPE',Line,iRc)
    call ChkIfKey()
  end if
  !---  Process NONE command
  if (Key('NONE')) then
    if (DBG) write(u6,*) ' Non-equilibrium response'
    NonEq = .true.
    call SetPos(LUInput,'NONE',Line,iRc)
    call ChkIfKey()
  end if
  !---  Process RFRO command
  if (Key('RFRO')) then
    call SetPos(LUInput,'RFRO',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    Line = Get_Ln(LUInput)
    Line(80:80) = '0'
    JPCMROOT = 0
    i_All = 0
    ReadStatus = ' Failure reading IPCMROOT after RFROOT keyword.'
    if (nToken(Line) >= 3) then
      read(Line,*,iostat=istatus) IPCMROOT,JPCMROOT,i_All
      if (istatus < 0) then
        call Error(2)
        return
      else if (istatus > 0) then
        call Error(3)
        return
      end if
      call DWSol_DWRO(LuInput,IPCMROOT,i_All)
    else
      read(Line,*,iostat=istatus) IPCMROOT
      if (istatus < 0) then
        call Error(2)
        return
      else if (istatus > 0) then
        call Error(3)
        return
      end if
    end if
    ReadStatus = ' O.K. reading IPCMROOT after RFROOT keyword.'

    ! Check that the root value is not changed explicitly by input
    jPCMRoot = iPCMRoot
    call Qpg_iScalar('RF0CASSCF root',Exists)
    if (Exists) then
      call Get_iScalar('RF0CASSCF root',jPCMRoot)
    else
      call Put_iScalar('RF0CASSCF root',iPCMRoot)
    end if

    if (jPCMRoot /= iPCMRoot) then
      !write(u6,*) 'iPCMRoot changed by explicitly by input.'

      ! Value changed explicitly by input. Accept the new value.

      call Put_iScalar('RF CASSCF root',iPCMRoot)
      call Put_iScalar('RF0CASSCF root',iPCMRoot)
      call Qpg_dArray('RF CASSCF Vector',Exists,mConf)
      if (Exists) then
        call mma_allocate(RF,mConf,Label='RF')
        RF(:) = Zero
        call Put_dArray('RF CASSCF Vector',RF,mConf)
        call mma_deallocate(RF)
      end if

    else

      ! Value not explicitly changed by input. See if it is
      ! changed by the RunFile, if it exists there.

      call Qpg_iScalar('RF CASSCF root',Exists)
      if (Exists) then
        call Get_iScalar('RF CASSCF root',iPCMRoot)
      else
        call Put_iScalar('RF CASSCF root',iPCMRoot)
      end if

    end if

    if (DBG) then
      write(u6,*) ' RFROOT command was given.'
      write(u6,*) ' Response field for root number ',IPCMROOT
    end if
    call ChkIfKey()
  end if
  !---  Process CIRF command
  if (Key('CIRF')) then
    call SetPos(LUInput,'CIRF',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' O.K. reading after keyword CIRF.'
    read(LUInput,*,iostat=istatus) ICIRFROOT
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. reading after keyword CIRF.'
    call ChkIfKey()
    if (DBG) then
      write(u6,*) ' CIRFROOT command was given.'
      write(u6,*) ' Response field will follow CISE root: ',ICIRFROOT
    end if
  end if
  !---------------------------------------------------------------------
  if (Key('MCM7')) then
#   ifndef _HDF5_
    call WarningMessage(2,'MCM7 is given in the input, please make sure to compile Molcas with HDF5 support.')
    call Error(4)
    return
#   endif
    MCM7 = .true.
    DoNECI = .true.  ! needed to initialise FCIQMC
    totalwalkers = 20000
    RDMsampling%start = 20000
    RDMsampling%n_samples = 10000
    RDMsampling%step = 100
    if (DBG) then
      write(u6,*) 'M7 CASSCF activated.'
      write(u6,*) 'Decoupled mode not implemented.'
      write(u6,*) 'Ignore automatically generated FciInp!'
    end if
  end if
  !---------------------------------------------------------------------
  if (Key('NDPT')) then
    tNonDiagStochPT2 = .true.
    IPT2 = 1     ! flag for SXCTL
    if (Key('SUPS')) then
      write(u6,*) 'SUPSymmetry incompatible with NDPT.'
      call Abend()
    end if
    if (Key('PPT2')) then
      write(u6,*) 'Non-diagonal PT2 incompatible with PPT2.'
      call Abend()
    end if
    if (DBG) then
      write(u6,*) 'stoch.-PT2 will be prepared in the current basis.'
      write(u6,*) 'Act. Space Fock matrix will be dumped.'
    end if
  end if
  !---------------------------------------------------------------------
  if (Key('PPT2')) then
    tPrepStochCASPT2 = .true.
    iOrbTyp = 2  ! pseudo-canonical orbitals
    IPT2 = 1     ! flag for SXCTL
    if (Key('SUPS')) then
      write(u6,*) 'SUPSymmetry incompatible with PPT2.'
      call Abend()
    end if
    if (Key('NDPT')) then
      write(u6,*) 'Non-diagonal PT2 incompatible with PPT2.'
      call Abend()
    end if
    if (DBG) then
      write(u6,*) 'Transforming final orbitals into pseudo-canonical.'
      write(u6,*) 'Act. Space Fock matrix will be dumped.'
    end if
  end if
  !---------------------------------------------------------------------
  if (Key('RGRA')) then
    tRootGrad = .true.
    if (DBG) write(u6,*) 'Orbital gradient for each root is printed.'
  end if
  !---  Process SSCR command
  if (Key('SSCR')) then
    if (DBG) write(u6,*) ' SSCR command was given.'
    call setpos(luinput,'SSCR',line,irc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    line = get_ln(luinput)
    line(80:80) = '0'
    ReadStatus = ' Failure reading after SSCR keyword.'
    read(line,*,iostat=istatus) norbs,same_orbs
    if (istatus /= 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K reading after SSCR keyword.'

    if (norbs >= mxOrb) then
      write(u6,'(a)',advance='no') 'SSCR error:'
      write(u6,*) 'number of spatial orbitals exceeds maximum'
      write(u6,'(a,i4)') 'norbs = ',norbs
      write(u6,'(a)') new_line('a')
      call abend()
    end if

    call mma_allocate(orb_range_p,norbs)
    call mma_allocate(orb_range_q,norbs)

    if (same_orbs /= 1) then
      Line = Get_Ln(LUInput)
      readstatus = ' failure reading after SSCR keyword.'
      read(Line,*) (orb_range_p(i),i=1,norbs)
      Line = Get_Ln(LUInput)
      read(Line,*) (orb_range_q(j),j=1,norbs)

      if (size(orb_range_p) /= size(orb_range_q)) then
        write(u6,'(a)',advance='no') 'SSCR error:'
        write(u6,*) 'numbers of spatial orbitals do not match'
        write(u6,*) 'orb_range_p has length ',size(orb_range_p)
        write(u6,*) 'orb_range_q has length ',size(orb_range_q)
        write(u6,'(a)') new_line('a')
        call abend()
      end if

      do i=1,norbs
        do j=1,norbs
          if (i < j) then
            if (orb_range_p(i) == orb_range_p(j)) then
              write(u6,'(a)',advance='no') 'SSCR error:'
              write(u6,*) 'first range contains duplicates.'
              write(u6,'(*(i4))') orb_range_p
              write(u6,'(a)') new_line('a')
              call abend()
            end if
            if (orb_range_q(i) == orb_range_q(j)) then
              write(u6,'(a)',advance='no') 'SSCR error:'
              write(u6,*) 'second range contains duplicates.'
              write(u6,'(*(i4))') orb_range_q
              write(u6,'(a)') new_line('a')
              call abend()
            end if
          end if
        end do
      end do
    else
      do i=1,norbs
      end do
      orb_range_p(1:norbs) = [(i,i=1,norbs)]
      orb_range_q(1:norbs) = [(i,i=1,norbs)]
    end if
    call ChkIfKey()
  end if
  !---  Process STAV command
  if (Key('STAV') .and. Key('CIRO')) then
    call WarningMessage(1,'STAVERAGE and CIROOT are incompatible.;The STAVERAGE command will be ignored.')
    call SetKey('STAV',.false.)
  end if
  if (Key('STAV')) then
    if (DBG) write(u6,*) ' STAVERAGE command was given.'
    call SetPos(LUInput,'STAV',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    Line = Get_Ln(LUInput)
    ReadStatus = ' Failure reading spin after STAVERAGE keyword.'
    read(Line,*,iostat=istatus) NROOTS
    if (istatus > 0) then
      call Error(3)
      return
    end if
    if (NROOTS > MXROOT) then
      write(u6,*) 'Error: number of roots exceeds maximum'
      write(u6,*) 'NROOTS = ',NROOTS
      write(u6,*) 'MXROOT = ',MXROOT
      call AbEnd()
    end if
    ReadStatus = ' O.K. reading spin after STAVERAGE keyword.'
    LROOTS = NROOTS
    iroot(1:NROOTS) = [(i,i=1,NROOTS)]
    WEIGHT(1:NROOTS) = One/real(NROOTS,kind=wp)
    if (DBG) then
      write(u6,*) ' Nr of roots in CI: LROOTS=',LROOTS
      write(u6,*) ' Nr of roots optimized by super-CI: NROOTS=',NROOTS
      write(u6,*) ' (Equal-weighted)'
    end if
    call ChkIfKey()
  end if
  !---  Process CIRO command
  if (DBG) write(u6,*) ' Check for CIROOTS command.'
  if (Key('CIRO')) then
    if (DBG) write(u6,*) ' CIROOTS command was given.'
    call SetPos(LUInput,'CIRO',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    Line = Get_Ln(LUInput)
    !BOR.. Modification 001011
    Line(80:80) = '0'
    ReadStatus = ' Failure reading after CIROOTS keyword.'
    read(Line,*,iostat=istatus) NROOTS,LROOTS,i_All
    if (istatus /= 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K reading after CIROOTS keyword.'
    if (NROOTS > MXROOT) then
      write(u6,*) 'Error: number of roots exceeds maximum'
      write(u6,*) 'NROOTS = ',NROOTS
      write(u6,*) 'MXROOT = ',MXROOT
      call AbEnd()
    end if
    if (i_All == 1) then
      iroot(1:NROOTS) = [(i,i=1,NROOTS)]
      WEIGHT(1:NROOTS) = One/real(NROOTS,kind=wp)
    else
      !BOR.. End modification 001011
      Line = Get_Ln(LUInput)
      ReadStatus = ' Failure reading after CIROOTS keyword.'
      read(Line,*,iostat=istatus) (IROOT(I),I=1,NROOTS)
    if (istatus /= 0) then
      call Error(3)
      return
    end if
      ReadStatus = ' O.K.after CIROOTS keyword.'
      WEIGHT(:) = Zero
      if (NROOTS == 1) then
        WEIGHT(1) = One
      else
        call mma_allocate(Temp1,NROOTS,Label='Temp1')
        Line = Get_Ln(LUInput)
        ReadStatus = ' Failure reading after CIROOTS keyword.'
        read(Line,*,iostat=istatus) (Temp1(i),i=1,NROOTS)
        if (istatus > 0) then
          call Error(3)
          return
        end if
        ReadStatus = ' O.K.after CIROOTS keyword.'
        iSum = sum(Temp1(1:NROOTS))
        WEIGHT(1:NROOTS) = real(Temp1(1:NROOTS),kind=wp)/real(iSum,kind=wp)
        call mma_deallocate(Temp1)
      end if
    end if
    if (DBG) then
      write(u6,*) ' Nr of roots in CI: LROOTS=',LROOTS
      write(u6,*) ' Nr of roots optimized by super-CI: NROOTS=',NROOTS
      if (i_All == 1) then
        write(u6,*) ' (Equal-weighted)'
      else
        write(u6,*) ' Weights:'
        do i1=1,NROOTS,10
          i2 = min(NROOTS,i1+9)
          write(u6,'(1x,10f8.4)') (WEIGHT(i),i=i1,i2)
        end do
      end if
    end if
    call ChkIfKey()
  end if
  !---  Process RLXR command
  if (Key('RLXR')) then
    call SetPos(LUInput,'RLXR',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' O.K. reading after keyword RLXR.'
    read(LUInput,*,iostat=istatus) IRLXROOT
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. reading after keyword RLXR.'
    call ChkIfKey()
    if (DBG) then
      write(u6,*) ' RLXROOT command was given.'
      write(u6,*) ' State for SLAPAF to handle: ',IRLXROOT
    end if
    if (.not. any(iRoot(1:LROOTS) == IRLXROOT)) then
      write(u6,*) ' The selected root is not among those available.'
      call AbEnd()
    end if
  end if

  !IgorS 29-4-2010 Begin
  !---  Process MDRL command
  if (Key('MDRL')) then
    if (Key('RLXR')) then
      write(u6,*) ' RLXROOT keyword was given before MDRLXROOT.'
      write(u6,*) ' Since these keywords are mutually exclusive'
      write(u6,*) ' please check the input and read the manual.'
      call Error(2)
      return
    end if
    call SetPos(LUInput,'MDRL',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    call Qpg_iScalar('Relax CASSCF root',RlxRCheck)
    if (RlxRCheck) then
      call Get_iScalar('Relax CASSCF root',iRlxRoot)
      if (DBG) then
        write(u6,*) ' An existing relax root was found.'
        write(u6,*) ' The MDRLXROOT value is ignored.'
      end if
    else
      ReadStatus = ' Failure reading relaxroot number IRLXROOT.'
      read(LUInput,*,iostat=istatus) IRLXROOT
      if (istatus < 0) then
        call Error(2)
        return
      else if (istatus > 0) then
        call Error(3)
        return
      end if
      ReadStatus = ' O.K. after reading relaxroot number IRLXROOT.'
      call ChkIfKey()
      if (.not. any(iRoot(1:LROOTS) == IRLXROOT)) then
        write(u6,*) ' The selected root is not among those available.'
        call AbEnd()
      end if
    end if
    if (DBG) then
      write(u6,*) ' MDRLxroot command was given.'
      write(u6,*) ' DYNAMIX will follow the root: ',IRLXROOT
    end if
  end if
  !IgorS End

  !---  Process CISE command
  if (DBG) write(u6,*) ' Check for CISELECT command.'
  if (Key('CISE')) then
    if (DBG) then
      write(u6,*) ' CISELECT keyword was given.'
      write(u6,*) ' This input is awkward. Let''s find up'
      write(u6,*) ' a better way to do things.'
    end if
    ICICH = 1
    call SetPos(LUInput,'CISE',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    do i=1,NROOTS
      ReadStatus = ' Failure reading after CISELECT keyword.'
      read(LUInput,*,iostat=istatus) kRef
      if (istatus < 0) then
        call Error(2)
        return
      else if (istatus > 0) then
        call Error(3)
        return
      end if
      ReadStatus = ' O.K. reading after CISELECT keyword.'
      if (kRef > mxRef) then
        call WarningMessage(1,'CISElect input is wrong.')
        write(u6,*) 'Number of CSF''s in CiSelect is out of bounds'
        write(u6,'(a,i3,a,i3)') 'Specified:',kRef,', Max is',mxRef
        write(u6,'(a,i3)') 'Standard fixup, value set to',mxRef
        kRef = mxRef
      end if
      ReadStatus = ' Failure reading after CISELECT keyword.'
      read(LUInput,*,iostat=istatus) (ICI(i,iRef),iRef=1,kRef)
      if (istatus < 0) then
        call Error(2)
        return
      else if (istatus > 0) then
        call Error(3)
        return
      end if
      read(LUInput,*,iostat=istatus) (CCI(i,iRef),iRef=1,kRef)
      if (istatus < 0) then
        call Error(2)
        return
      else if (istatus > 0) then
        call Error(3)
        return
      end if
      ReadStatus = ' O.K. reading after CISELECT keyword.'
      dSum = sum(CCI(i,1:kRef)**2)
      CCI(i,1:kRef) = CCI(i,1:kRef)/sqrt(dSum)
      CCI(i,kRef+1:mxRef) = Zero
      ICI(i,kRef+1:mxRef) = 0
    end do
    call ChkIfKey()
  end if

  !---  Process ALPH command
  if (Key('ALPH')) then
    if (DBG) write(u6,*) ' The ALPH keyword was used.'
    call SetPos(LUInput,'ALPH',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading data after ALPH keyword.'
    read(LUInput,*,iostat=istatus) iAlphaBeta
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading data after ALPH keyword.'
    if (iAlphaBeta < 0) iAlphaBeta = -1
    if (iAlphaBeta > 0) iAlphaBeta = 1
    if (DBG) then
      if (iAlphaBeta == 1) write(u6,*) ' Read alpha orbitals from UHF'
      if (iAlphaBeta == -1) write(u6,*) ' Read beta orbitals from UHF'
    end if
    call ChkIfKey()
  end if

  ! =========   Input source for orbitals: =============================
  ! INVEC=0 is used to indicate if any source of orbitals has been
  ! identified.
  ! iOrbData=0 is used to indicate if any source of orbital type
  ! information -- inactive, ras1, etc -- has been identified.

  if (DBG) then
    write(u6,*) ' Check for input source of orbitals.'
    write(u6,*) 'Key("CORE"),Key("JOBI"),Key("LUMO"):',Key('CORE'),Key('JOBI'),Key('LUMO')
  end if
  ! Handle multiple keywords:
  if (Key('CORE')) then
    call SetKey('JOBI',.false.)
    call SetKey('LUMO',.false.)
  else if (Key('JOBI')) then
    call SetKey('LUMO',.false.)
  end if
  ! Only one of these should have been selected.
  if (DBG) write(u6,*) 'Key("CORE"),Key("JOBI"),Key("LUMO"):',Key('CORE'),Key('JOBI'),Key('LUMO')

  ! CORE is probably becoming obsolete.
  !---  Process CORE command
  if (Key('CORE')) then
    if (DBG) write(u6,*) ' CORE command was used.'
    if (IPRLEV >= VERBOSE) write(u6,*) ' Start orbitals will be computed from scratch.'
    INVEC = 1
  end if

  !---  Process JOBI command
  if (Key('JOBI')) then
    if (DBG) write(u6,*) ' JOBIPH command was used.'
    call f_Inquire('JOBOLD',lExists)
    if (lExists) then
      if (IPRLEV >= VERBOSE) write(u6,*) ' Orbitals will be taken from the old jobiph file.'
      INVEC = 3
    else
      call f_Inquire(IPHNAME,lExists)
      if (lExists) then
        if (IPRLEV >= VERBOSE) write(u6,*) ' Orbitals will be taken from the jobiph file.'
        INVEC = 3
      end if
    end if
    if (INVEC == 0) then
      if (IPRLEV >= TERSE) then
        call WarningMessage(2,'JOBIPH input is wrong.')
        write(u6,*) ' Keyword JOBIPH was used, but the ''JOBOLD'' file'
        write(u6,*) ' does not exist. The ''JOBIPH'' file named'
        write(u6,*) IPHNAME
        write(u6,*) 'also does not exist. This is a fatal error.'
        call Error(4)
        return
      end if
    end if
  end if

  !---  Process H5OR command
  if (Key('H5OR')) then
#   ifdef _HDF5_
    call SetKey('LUMO',.false.)
    call SetKey('TYPE',.false.)
    iOverwr = merge(1,0,any([Key('RAS1'),Key('RAS2'),Key('RAS3'),Key('FROZ'),Key('INAC'),Key('DELE')]))
    mh5id = mh5_open_file_r(StartOrbFile)
    ! read basic attributes
    call mh5_fetch_attr(mh5id,'NSYM',NSYM_L)
    if (nsym /= nsym_l) then
      write(u6,*) 'Number of symmetries on HDF5 file does not'
      write(u6,*) 'match the number of symmetries on the'
      write(u6,*) 'RunFile, calculation will stop now.'
      call Quit(_RC_INPUT_ERROR_)
    end if
    call mh5_fetch_attr(mh5id,'NBAS',NBAS_L)
    ierr = 0
    do isym=1,nsym
      if (nbas(isym) /= nbas_l(isym)) ierr = 1
    end do
    if (ierr == 1) then
      write(u6,*) 'Number of basis functions on HDF5 file does not'
      write(u6,*) 'match the number of basis functions on the'
      write(u6,*) 'RunFile, calculation will stop now.'
      call Quit(_RC_INPUT_ERROR_)
    end if
    ! orbitals available?
    select case (iAlphaBeta)
      case (1)
        Line = 'MO_ALPHA_VECTORS'
      case (-1)
        Line = 'MO_BETA_VECTORS'
      case default
        Line = 'MO_VECTORS'
    end select
    if (mh5_exists_dset(mh5id,trim(Line))) inVec = 4
    ! typeindex data available?
    select case (iAlphaBeta)
      case (1)
        Line = 'MO_ALPHA_TYPEINDICES'
      case (-1)
        Line = 'MO_BETA_TYPEINDICES'
      case default
        Line = 'MO_TYPEINDICES'
    end select
    if (mh5_exists_dset(mh5id,trim(Line))) then
      iOrbData = 3
      call mma_allocate(typestring,sum(nbas(1:nsym)))
      call mh5_fetch_dset(mh5id,trim(Line),typestring)
      call tpstr2orb(nsym_l,nbas_l,typestring,nfro_l,nish_l,nrs1_l,nrs2_l,nrs3_l,nssh_l,ndel_l)
      call mma_deallocate(typestring)
    end if
    ! CI requested?
    if (Key('H5CI')) then
      if (mh5_exists_dset(mh5id,'CI_VECTORS') .and. mh5_exists_attr(mh5id,'NROOTS')) then
        call mh5_fetch_attr(mh5id,'NROOTS',lroots_l)
        if (lroots_l /= lroots) then
          write(u6,*) 'Number of CI roots on file does not'
          write(u6,*) 'match the number requested by the user,'
          write(u6,*) 'so no CI vectors will be read from'
          write(u6,*) StartOrbFile
          iCIRST = 0
        else
          iCIRST = 1
        end if
      else
        write(u6,*) 'The required fields CI_VECTORS and/or NROOTS'
        write(u6,*) 'are missing from the HDF5 file supplied by'
        write(u6,*) 'the user. As a result, to continue,'
        write(u6,*) 'no CI vectors will be read from'
        write(u6,*) StartOrbFile
        iCIRST = 0
      end if
    end if
    call mh5_close_file(mh5id)
#   else
    write(u6,*) 'The format of the start orbital file was'
    write(u6,*) 'specified by the user as HDF5, but this'
    write(u6,*) 'is not implemented in this installation.'
    call Quit(_RC_INPUT_ERROR_)
#   endif
  end if

  !---  Process LUMO command
  InfoLbl = '        '
  if (Key('LUMO')) then
    if (DBG) then
      write(u6,*) ' LUMORB command was used.'
      write(u6,*) ' Name of orbital file, StartOrbFile='
      write(u6,*) StartOrbFile
    end if
    call ChkVec(StartOrbFile,inporb_version,NSYM_L,NBAS_L,NORB_L,InfoLbl,iRc1)
    if (iRc1 /= _RC_ALL_IS_WELL_) then
      if (IPRLEV >= TERSE) then
        call WarningMessage(1,'LUMORB input error.')
        write(u6,*) ' Keyword LUMORB used with file name StartOrbFile='
        write(u6,*) StartOrbFile
        write(u6,*) ' but that file cannot be used. Perhaps it does'
        write(u6,*) ' not exist?'
      end if
      iOverWr = 1
    else
      if (DBG) then
        write(u6,*) ' The file may be used to read input orbitals.'
        write(u6,*) ' It is of type INPORB version ',inporb_version
        write(u6,*) ' The information label is: ',InfoLbl
      end if
      ! Check that symmetry and basis function information match the runfile:
      IERR = 0
      if (NSYM /= NSYM_L) then
        IERR = 1
      else
        do ISYM=1,NSYM
          if (NBAS(ISYM) /= NBAS_L(ISYM)) IERR = 2
        end do
      end if
      if (IERR /= 0) then
        call WarningMessage(2,'Mismatch in start orbital file.')
        write(u6,*) ' ERROR: Start orbital file name is '
        write(u6,*) StartOrbFile
        write(u6,*) ' That file is a valid orbital file.'
        write(u6,*) ' Version:',inporb_version
        write(u6,*) ' But some information does not match.'
        if (IERR == 1) then
          write(u6,*) ' In the file, nr of symmetries is =',NSYM_L
          write(u6,*) ' but according to the runfile, it is=',NSYM
        else if (IERR == 2) then
          write(u6,*) ' In the file, nr of basis functions/symm is'
          write(u6,'(1x,8I5)') (NBAS_L(I),I=1,NSYM)
          write(u6,*) ' but RUNFILE claims it is'
          write(u6,'(1x,8I5)') (NBAS(I),I=1,NSYM)
        end if
        write(u6,*) ' Is it an old file left in workspace by mistake?'
        call Error(4)
        return
      end if

      ! This also implies that information on orbital types could be
      ! taken from typeindex on orbital file:
      if (('I' .in. to_upper(trim(InfoLbl))) .and. &
          (.not. any([Key('RAS1'),Key('RAS2'),Key('RAS3'),Key('FROZ'),Key('INAC'),Key('DELE')]))) then
        iOrbData = 3
        iOverWr = 0
        if (DBG) then
          write(u6,*) ' This means we may take orbital specifications'
          write(u6,*) ' from the file, so set iOrbData=3, iOverWr=0.'
          write(u6,*) ' The orbital spaces are read from typeindex.'
        end if
        ! We will also take the opportunity to find the orbital spaces size
        ! according to typeindex, for possible need below:
        call mma_allocate(iType,mxOrb,Label='iType')
        LuStartOrb = 19
        call RdVec(StartOrbFile,LuStartOrb,'IA',NSYM_L,NBAS_L,NBAS_L,Dummy,Dummy,Dummy,iType,myTitle,0,iErr)
        call tpidx2orb(NSYM_L,NBAS_L,iType,NFRO_L,NISH_L,NRS1_L,NRS2_L,NRS3_L,NSSH_L,NDEL_L)
        call mma_deallocate(iType)
        if (DBG) then
          write(u6,*) ' From RDTPIDX, we get:'
          write(u6,'(1x,A16,8I4)') ' NBAS_L:',(NBAS_L(I),I=1,NSYM_L)
          write(u6,'(1x,A16,8I4)') ' NORB_L:',(NORB_L(I),I=1,NSYM_L)
          write(u6,'(1x,A16,8I4)') ' NFRO_L:',(NFRO_L(I),I=1,NSYM_L)
          write(u6,'(1x,A16,8I4)') ' NISH_L:',(NISH_L(I),I=1,NSYM_L)
          write(u6,'(1x,A16,8I4)') ' NRS1_L:',(NRS1_L(I),I=1,NSYM_L)
          write(u6,'(1x,A16,8I4)') ' NRS2_L:',(NRS2_L(I),I=1,NSYM_L)
          write(u6,'(1x,A16,8I4)') ' NRS3_L:',(NRS3_L(I),I=1,NSYM_L)
          write(u6,'(1x,A16,8I4)') ' NSSH_L:',(NSSH_L(I),I=1,NSYM_L)
          write(u6,'(1x,A16,8I4)') ' NDEL_L:',(NDEL_L(I),I=1,NSYM_L)
        end if
      else
        iOverWr = 1
      end if
    end if
    INVEC = 2
  end if
  if (DBG) then
    write(u6,*) ' The INVEC    code is now',INVEC
    write(u6,*) ' The iOrbData code is now',iOrbData
  end if

  ! ====================================================================

  !---  Process TYPEindex command
  if (DBG) write(u6,*) ' Was TYPEINDEX requested?'
  if (Key('TYPE')) then
    call SetPos(LUInput,'TYPE',Line,iRc)
    call ChkIfKey()
    if (DBG) then
      write(u6,*) ' TYPEINDEX command was used.'
      write(u6,*) ' The size of orbital spaces should be read from'
      write(u6,*) ' typeindex in starting orbital file.'
    end if
    if ((index(InfoLbl,'i') > 0) .or. (index(InfoLbl,'I') > 0)) then
      iOrbData = 3
      iOverwr = 0
      if (IPRLEV >= VERBOSE) write(u6,*) ' Orbital specification will be taken from orbital file'
      call mma_allocate(iType,mxOrb,Label='iType')
      LuStartOrb = 19
      call RdVec(StartOrbFile,LuStartOrb,'IA',NSYM_L,NBAS_L,NBAS_L,Dummy,Dummy,Dummy,iType,myTitle,0,iErr)
      call tpidx2orb(NSYM_L,NBAS_L,iType,NFRO_L,NISH_L,NRS1_L,NRS2_L,NRS3_L,NSSH_L,NDEL_L)
      call mma_deallocate(iType)
      IERR = 0
      if (NSYM_L /= NSYM) IERR = 1
      if (IERR == 0) then
        do ISYM=1,NSYM
          if (NBAS(ISYM) /= NBAS_L(ISYM)) IERR = 1
        end do
      end if
      if ((IERR /= 0) .or. DBG) then
        write(u6,*)
        write(u6,'(6X,A)') 'Specifications read from runfile:'
        write(u6,'(6X,A)') '----------------------------------------'
        write(u6,*)
        write(u6,'(6X,A,T47,8I4)') 'Symmetry species',(iSym,iSym=1,NSYM)
        write(u6,'(6X,A,T47,8I4)') 'Number of basis functions',(NBAS(iSym),iSym=1,NSYM)
        write(u6,*)
        write(u6,'(6X,A)') 'Specifications read from orbital file:'
        write(u6,'(6X,A)') '----------------------------------------'
        write(u6,*)
        write(u6,'(6X,A,T47,8I4)') 'Symmetry species',(iSym,iSym=1,NSYM_L)
        write(u6,'(6X,A,T47,8I4)') 'Frozen orbitals',(NFRO_L(iSym),iSym=1,NSYM_L)
        write(u6,'(6X,A,T47,8I4)') 'Inactive orbitals',(NISH_L(iSym),iSym=1,NSYM_L)
        write(u6,'(6X,A,T47,8I4)') 'RAS1 orbitals',(NRS1_L(iSym),iSym=1,NSYM_L)
        write(u6,'(6X,A,T47,8I4)') 'RAS2 orbitals',(NRS2_L(iSym),iSym=1,NSYM_L)
        write(u6,'(6X,A,T47,8I4)') 'RAS3 orbitals',(NRS3_L(iSym),iSym=1,NSYM_L)
        write(u6,'(6X,A,T47,8I4)') 'Secondary orbitals',(NSSH_L(iSym),iSym=1,NSYM_L)
        write(u6,'(6X,A,T47,8I4)') 'Deleted orbitals',(NDEL_L(iSym),iSym=1,NSYM_L)
        write(u6,'(6X,A,T47,8I4)') 'Number of basis functions',(NBAS_L(iSym),iSym=1,NSYM_L)
        write(u6,*)
      end if
      if ((IERR /= 0) .and. (IPRLEV >= TERSE)) then
        write(u6,*) ' Orbital specifications were to be read from'
        write(u6,*) ' orbital file, but there is mismatch with'
        write(u6,*) ' some data on the runfile!'
        write(u6,*) ' Orbital file name is:'
        write(u6,*) StartOrbFile
        iOrbData = 0
      else
        NFRO(1:NSYM) = NFRO_L(1:NSYM)
        NISH(1:NSYM) = NISH_L(1:NSYM)
        NRS1(1:NSYM) = NRS1_L(1:NSYM)
        NRS2(1:NSYM) = NRS2_L(1:NSYM)
        NRS3(1:NSYM) = NRS3_L(1:NSYM)
        NSSH(1:NSYM) = NSSH_L(1:NSYM)
        NDEL(1:NSYM) = NDEL_L(1:NSYM)
      end if
    else
      if (DBG) write(u6,*) ' Typeindex cannot be read!'
    end if
  end if
  if (DBG) write(u6,*) ' The iOrbData code is now',iOrbData

  ! ====================================================================
  ! Explicit orbital sizes input:
  ! Save a copy of current iorbdata first:
  iod_save = iorbdata
  !---  Process FROZ command
  if (Key('FROZ')) then
    if (DBG) write(u6,*) ' FROZEN keyword was given.'
    call SetPos(LUInput,'FROZ',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading after FROZEN keyword.'
    read(LUInput,*,iostat=istatus) (NFRO(ISYM),ISYM=1,NSYM)
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. reading after FROZEN keyword.'
    call Get_iScalar('nSym',i)
    call Put_iArray('nFro',nFro,i)
    if (DBG) then
      write(u6,*) ' Nr of Frozen orbitals requested:'
      write(u6,'(1x,8i5)') (NFRO(i),i=1,NSYM)
    end if
    IORBDATA = 1
  end if
  !---  Process INAC command
  if (Key('INAC')) then
    if (DBG) write(u6,*) ' INACTIVE keyword was given.'
    call SetPos(LUInput,'INAC',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading after INACTIVE keyword.'
    read(LUInput,*,iostat=istatus) (NISH(ISYM),ISYM=1,NSYM)
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. reading after INACTIVE keyword.'
    if (DBG) then
      write(u6,*) ' Nr of Inactive orbitals requested:'
      write(u6,'(1x,8i5)') (NISH(i),i=1,NSYM)
    end if
    IORBDATA = 1
  end if

  !---  Process RAS1 command
  if (Key('RAS1')) then
    if (DBG) write(u6,*) ' RAS1 keyword was given.'
    call SetPos(LUInput,'RAS1',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading after RAS1 keyword.'
    read(LUInput,*,iostat=istatus) (NRS1(ISYM),ISYM=1,NSYM)
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. reading after RAS1 keyword.'
    if (DBG) then
      write(u6,*) ' Nr of RAS1 orbitals requested:'
      write(u6,'(1x,8i5)') (NRS1(i),i=1,NSYM)
    end if
    IORBDATA = 1
  end if

  !---  Process RAS2 command
  if (Key('RAS2')) then
    if (DBG) write(u6,*) ' RAS2 keyword was given.'
    call SetPos(LUInput,'RAS2',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading after RAS2 keyword.'
    read(LUInput,*,iostat=istatus) (NRS2(ISYM),ISYM=1,NSYM)
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. reading after RAS2 keyword.'
    if (DBG) then
      write(u6,*) ' Nr of Ras2 orbitals requested:'
      write(u6,'(1x,8i5)') (NRS2(i),i=1,NSYM)
    end if
    IORBDATA = 1
  end if

  !---  Process RAS3 command
  if (Key('RAS3')) then
    if (DBG) write(u6,*) ' RAS3 keyword was given.'
    call SetPos(LUInput,'RAS3',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading after RAS3 keyword.'
    read(LUInput,*,iostat=istatus) (NRS3(ISYM),ISYM=1,NSYM)
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. reading after RAS3 keyword.'
    if (DBG) then
      write(u6,*) ' Nr of Ras3 orbitals requested:'
      write(u6,'(1x,8i5)') (NRS3(i),i=1,NSYM)
    end if
    IORBDATA = 1
  end if

  !---  Process GASS command
  if (Key('GASS')) then
    if (DBG) write(u6,*) 'GAS is actived'
    call setpos(luinput,'GASS',line,irc)
    if (irc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    read(luinput,*,iostat=istatus) NGAS
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    do igas=1,ngas
      read(luinput,*,iostat=istatus) (ngssh(igas,isym),isym=1,nsym)
      if (istatus < 0) then
        call Error(2)
        return
      else if (istatus > 0) then
        call Error(3)
        return
      end if
      read(luinput,*,iostat=istatus) (igsoccx(igas,mm),mm=1,2)
      if (istatus < 0) then
        call Error(2)
        return
      else if (istatus > 0) then
        call Error(3)
        return
      end if
    end do
    iDoGas = .true.
    iorbdata = 1
  end if

  !---  Process DELE command
  if (Key('DELE')) then
    if (DBG) write(u6,*) ' DELETED keyword was given.'
    call SetPos(LUInput,'DELE',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading after DELETED keyword.'
    read(LUInput,*,iostat=istatus) (NDEL(ISYM),ISYM=1,NSYM)
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading after DELETED keyword.'
    if (DBG) then
      write(u6,*) ' Nr of Deleted orbitals requested:'
      write(u6,'(1x,8i5)') (NDEL(i),i=1,NSYM)
    end if
    IORBDATA = 1
  end if

  if ((IORBDATA == 1) .and. (IPRLEV >= VERBOSE)) write(u6,*) ' Orbital specification was read from input.'
  if ((IOD_SAVE == 3) .and. (IORBDATA == 1)) then
    ! See if the input matches the values on file
    IERR = 0
    if (any(NFRO(1:NSYM) /= NFRO_L(1:NSYM))) IERR = 1
    if (any(NISH(1:NSYM) /= NISH_L(1:NSYM))) IERR = 1
    if (any(NRS1(1:NSYM) /= NRS1_L(1:NSYM))) IERR = 1
    if (any(NRS2(1:NSYM) /= NRS2_L(1:NSYM))) IERR = 1
    if (any(NRS3(1:NSYM) /= NRS3_L(1:NSYM))) IERR = 1
    if (any(NSSH(1:NSYM) /= NSSH_L(1:NSYM))) IERR = 1
    if (any(NDEL(1:NSYM) /= NDEL_L(1:NSYM))) IERR = 1
    if (IERR == 0) then
      if ((IORBDATA == 1) .and. (IPRLEV >= VERBOSE)) then
        write(u6,*) ' However, input matches the typeindex on the'
        write(u6,*) ' starting orbitals file. Therefore, accept'
        write(u6,*) ' the typeindex information for sorting.'
        iOrbData = 3
        iOverWr = 0
      end if
    end if
  end if
  if (IORBDATA == 3) then
    NFRO(1:NSYM) = NFRO_L(1:NSYM)
    NISH(1:NSYM) = NISH_L(1:NSYM)
    NRS1(1:NSYM) = NRS1_L(1:NSYM)
    NRS2(1:NSYM) = NRS2_L(1:NSYM)
    NRS3(1:NSYM) = NRS3_L(1:NSYM)
    NSSH(1:NSYM) = NSSH_L(1:NSYM)
    NDEL(1:NSYM) = NDEL_L(1:NSYM)
  end if
  ! ====================================================================
  ! If IORBDATA is still 0, lets hope there is information on the runfile.
  ! Exception: If this is a CIRESTART, it must be taken from the JOBIPH
  ! (or JOBOLD) file.
  if (IORBDATA == 0) then
    if (IPRLEV >= VERBOSE) write(u6,*) ' No explicit orbital specs in user input.'
    if (Key('CIRE')) then
      if (IPRLEV >= VERBOSE) then
        write(u6,*) ' This is a CIRESTART case, so take them from'
        write(u6,*) ' the JOBIPH or JOBOLD file.'
      end if
      IORBDATA = 2
      if (IPRLEV >= VERBOSE) write(u6,*) ' Orbital specs taken from JOBIPH or JOBOLD.'
      IAD19 = 0
      iJOB = 0
      call f_Inquire('JOBOLD',lExists)
      if (lExists) iJOB = 1
      if (JOBOLD <= 0) JOBOLD = 20
      if (iJOB == 1) then
        call DaName(JOBOLD,'JOBOLD')
      else
        if (IPRLEV >= TERSE) call WarningMessage(1,'JOBOLD not found, using JOBIPH.')
        if (JOBIPH > 0) then
          JOBOLD = JOBIPH
        else
          call DaName(JOBOLD,'JOBIPH')
        end if
      end if
      call IDaFile(JOBOLD,2,IADR19,10,IAD19)
      ! PAM Jan 2014 -- do not take POTNUC from JOBIPH; take it directly
      ! from runfile, where it was stored by seward.
      iAd19 = iAdr19(1)
      call WR_RASSCF_Info(JobOld,2,iAd19,NACTEL,ISPIN,NSYM,STSYM,NFRO,NISH,NASH,NDEL,NBAS,mxSym,lJobH1,(LenIn+8)*mxOrb,NCONF, &
                          lJobH2,2*72,JobTit,4*18*mxTit,POTNUCDUMMY,LROOTS,NROOTS,IROOT,mxRoot,NRS1,NRS2,NRS3,NHOLE1,NELEC3,IPT2, &
                          WEIGHT)
    else
      if (IPRLEV >= VERBOSE) then
        write(u6,*) ' This is not a CIRESTART case, so take them from'
        write(u6,*) ' the RUNFILE.'
      end if
      call qpg_iArray('nAsh',lExists,nItems)
      if (lExists .and. (nItems == nSym)) then
        call Get_iArray('nFro',nFro,nSym)
        call Get_iArray('nISh',nISh,nSym)
        call Get_iArray('nASh',nRS2,nSym)
        call Get_iArray('nDel',nDel,nSym)
        IORBDATA = 2
      end if
    end if
  end if
  ! complete orbital specifications
  do iSym=1,nSym
    if (.not. iDoGas) then
      nash(isym) = nrs1(isym)+nrs2(isym)+nrs3(isym)
    else
      NASH(ISYM) = sum(NGSSH(1:NGAS,ISYM))
    end if
    NORB(ISYM) = NBAS(ISYM)-NFRO(ISYM)-NDEL(ISYM)
    NSSH(ISYM) = NORB(ISYM)-NISH(ISYM)-NASH(ISYM)
  end do
  ! Related data for sizes, etc.
  NTOT = sum(NBAS(1:NSYM))
  NTOT1 = sum(nTri_Elem(NBAS(1:NSYM)))
  NTOT2 = sum(NBAS(1:NSYM)**2)
  NO2M = maxval(NBAS(1:NSYM)**2)
  NISHT = sum(NISH(1:NSYM))
  NASHT = sum(NASH(1:NSYM))
  NDELT = sum(NDEL(1:NSYM))
  NFROT = sum(NFRO(1:NSYM))
  NSEC = sum(NSSH(1:NSYM))
  NORBT = sum(NORB(1:NSYM))
  NTOT3 = sum(nTri_Elem(NORB(1:NSYM)))
  NTOTSP = sum(nTri_Elem(NASH(1:NSYM)))
  NTOT4 = sum(NORB(1:NSYM)**2)
  NRS1T = sum(NRS1(1:NSYM)) ! for RASSCF
  NRS2T = sum(NRS2(1:NSYM))
  NRS3T = sum(NRS3(1:NSYM))
  !NGSSH_tot(:) = Zero
  !do igas=1,ngas
  !  NGSSH_tot(igas) = sum(NGSSH(igas,1:NSYM))
  !end do
  !if (dbg) then
  !  write(u6,*) 'NGSSH_tot(igas):'
  !  write(u6,*) (NGSSH_tot(igas),igas=1,ngas)
  !end if
  NACPAR = nTri_Elem(NASHT)
  NACPR2 = nTri_Elem(NACPAR)
  ! NASHT is called NAC in some places:
  NAC = NASHT
  ! Same, NISHT, NIN:
  NIN = NISHT
  NFR = NFROT
  if (DBG) write(u6,*) ' The iOrbData code is now',iOrbData
  ! ====================================================================
  ! Compute effective nuclear charge.
  ! Identical to nr of protons for conventional basis sets only, not ECP.
  call Get_iScalar('Unique atoms',nNuc)
  call mma_allocate(ENC,nNuc,Label='ENC')
  call Get_dArray('Effective nuclear Charge',ENC,nNuc)
  TEffNChrg = Zero
  call mma_allocate(Stab,nNuc,Label='Stab')
  call Get_iArray('nStab',Stab,nNuc)
  do i=1,nNuc
    TEffNChrg = TEffNChrg+ENC(i)*real(nSym/Stab(i),kind=wp)
  end do
  call mma_deallocate(Stab)
  call mma_deallocate(ENC)
  if (DBG) write(u6,*) ' Effective nuclear charge is TEffNChrg=',TEffNChrg
  TotChrg = Zero
  if (DBG) write(u6,*) ' Set TotChrg=',TotChrg
  !---  Process NACT command
  if (Key('NACT')) then
    ! Cannot set the number of inactive orbitals if explicitly given or
    ! or if there is symmetry
    if (Key('CHAR') .and. (Key('INAC') .or. (NSYM > 1))) then
      if (IPRLEV >= TERSE) call WarningMessage(1,'CHARGE and NACTEL are only compatible if INACTIVE is not given and the '// &
                                               'symmetry group is C1. Hence the CHARGE command will be ignored.')
      call SetKey('CHAR',.false.)
    end if
    if (DBG) write(u6,*) ' NACTEL keyword was given.'
    call SetPos(LUInput,'NACT',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    Line = Get_Ln(LUInput)
    ReadStatus = ' Failure reading data after NACTEL keyword.'
    read(Line,*,iostat=istatus) NACTEL,NHOLE1,NELEC3
    if (istatus < 0) then
      read(Line,*,iostat=istatus) NACTEL
      if (istatus /= 0) then
        call Error(3)
        return
      end if
      NHOLE1 = 0
      NELEC3 = 0
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading data after NACTEL keyword.'
    if (DBG) write(u6,*) ' NACTEL,NHOLE1,NELEC3:',NACTEL,NHOLE1,NELEC3
    ! Only set total charge here if not explicitly given
    if (.not. Key('CHAR')) then
      TotChrg = TEffNChrg-real(Two*(NISHT+NFROT)+NACTEL)
      if (DBG) write(u6,*) ' TotChrg=',TotChrg
    end if
    call Put_iScalar('nActel',NACTEL)
    call ChkIfKey()
  end if
  !---  Process CHARGE command
  if (Key('CHAR')) then
    if (DBG) write(u6,*) ' CHARGE command was given.'
    call SetPos(LUInput,'CHAR',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    Line = Get_Ln(LUInput)
    ReadStatus = ' Failure reading charge after CHARGE keyword.'
    read(Line,*,iostat=istatus) NCHRG
    if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. reading charge after CHARGE keyword.'
    if (DBG) write(u6,*) ' Total charge is ',NCHRG
    TotChrg = real(NCHRG,kind=wp)
    call ChkIfKey()
    ! If both CHAR and NACT where given, set the inactive and secondary
    ! orbitals accordingly
    if (Key('NACT')) then
      NISHT_old = NISHT
      NISHT = (int(TEffNChrg-TotChrg+Half)-NACTEL)/2-NFROT
      NISH(1) = NISHT
      NIN = NISHT
      NSEC = NSEC-(NISHT-NISHT_old)
      NSSH(1) = NSEC
    end if
  end if
  ! The NINT function is unreliable on Cygwin gfortran, use INT:
  ! Nr of electrons should be positive integer, so this is probably safe:
  NACTEL = int(TEffNChrg-TotChrg+Half)-2*(NISHT+NFROT)
  call Put_iScalar('nActel',NACTEL)
  if (DBG) then
    write(u6,*) ' Compute NActEl from  other data:'
    write(u6,*) '     TEffNChrg=',TEffNChrg
    write(u6,*) '       TotChrg=',TotChrg
    write(u6,*) '       NISHT  =',NISHT
    write(u6,*) '       NFROT  =',NFROT
    write(u6,*) ' Resulting NActEl=',NActEl
  end if
  !---  Process RASSCF command
  if (Key('RASS')) then
    if (DBG) write(u6,*) ' RASSCF keyword was given.'
    call SetPos(LUInput,'RASS',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    Line = Get_Ln(LUInput)
    ReadStatus = ' Failure reading data after RASSCF keyword.'
    read(Line,*,iostat=istatus) NHOLE1,NELEC3
    if (istatus /= 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. reading data after RASSCF keyword.'
    call ChkIfKey()
  end if
  ! ====================================================================
  !---  Process SPIN command
  if (DBG) write(u6,*) ' Determine spin value:'
  if (Key('SPIN')) then
    if (DBG) write(u6,*) ' SPIN command was given.'
    call SetPos(LUInput,'SPIN',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    Line = Get_Ln(LUInput)
    ReadStatus = ' Failure reading spin after SPIN keyword.'
    read(Line,*,iostat=istatus) ISPIN
    if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. reading spin after SPIN keyword.'
    if (DBG) then
      write(u6,*) ' Spin multiplicity is ',ISPIN
      S = Half*real(ISPIN-1,kind=wp)
      write(u6,*) '  i.e. SPIN=',S
    end if
    call ChkIfKey()
  else
    ! If ISPIN has not been set, use some value that may have been set
    ! in runfile (e.g., UHF, or previous calculation)
    if (DBG) write(u6,*) ' No SPIN command was given.'
    call qpg_iscalar('ISPIN',lExists)
    if (lExists) then
      call get_iscalar('ISPIN',ISPIN)
      if (DBG) write(u6,*) ' Runfile has ISPIN=',ISPIN
    else
      if (DBG) write(u6,*) ' Runfile does not know ISPIN.'
      call qpg_dscalar('UHFSPIN',lExists)
      if (lExists) then
        call get_dscalar('UHFSPIN',SUHF)
        ISPIN = nint(One+Two*SUHF)
        if (DBG) write(u6,*) ' Runfile has UHFSPIN=',SUHF
      else
        ! or, last chance fallback, guess on singlet or doublet:
        if (DBG) write(u6,*) ' Runfile does not know UHFSPIN.'
        ISPIN = 1
        if (mod(NACTEL,2) /= 0) ISPIN = 2
        if (IPRLEV >= TERSE) then
          call WarningMessage(1,'Had to guess the spin.')
          write(u6,*) ' Warning: no input and no reliable source'
          write(u6,*) ' for the spin multiplicity.'
          write(u6,*) ' Guess ISPIN=',ISPIN
        end if
      end if
    end if
  end if
  call put_iscalar('ISPIN',ISPIN)
  ! If spin is zero, do not compute and print spin density:
  if (ISPIN == 1) ISPDEN = 0
  ! ====================================================================
  if (Key('DMPO')) DumpOnly = .true.
  ! ====================================================================
  if (Key('REOR')) then
    if (DBG) write(u6,*) 'Orbital Reordering (REOR) is activated'
    call setpos(luinput,'REOR',line,irc)
    if (irc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading ReOrFlag after REOR keyword.'
    read(luinput,*,iostat=istatus) ReOrFlag
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. reading ReOrFlag after REOR keyword.'
    if ((ReOrFlag < -1) .or. (ReOrFlag == 1)) then
      call WarningMessage(2,"Invalid flag for reordering. n==0: Don't reorder. n>=2: User defined permutation with n changed "// &
                          'elements. n==-1: Use GAS sorting scheme. ')
      call Error(4)
      return
    else if (ReOrFlag >= 2) then
      call mma_allocate(ReOrInp,ReOrFlag)
      ReadStatus = ' Failure reading ReOrInp after REOR keyword.'
      read(luinput,*,iostat=istatus) (ReOrInp(i),i=1,ReOrFlag)
      if (istatus < 0) then
        call Error(2)
        return
      else if (istatus > 0) then
        call Error(3)
        return
      end if
      ReadStatus = ' O.K. reading ReOrInp after REOR keyword.'
    else if ((ReOrFlag == -1) .and. (.not. Key('GASS'))) then
      call WarningMessage(2,'If GAS is not used, a permutation for orbital reordering has to be specified.')
      call Error(4)
      return
    end if
  end if
  if (Key('ORTH')) then
    call setpos(luinput,'ORTH',line,irc)
    if (irc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    read(luinput,*,iostat=istatus) ON_scheme_inp
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    uppercased = to_upper(trim(ON_scheme_inp))
    if ('CANO' .in. uppercased) then
      ON_scheme%val = ON_scheme_values%Canonical
    else if ('LOWD' .in. uppercased) then
      ON_scheme%val = ON_scheme_values%Lowdin
    else if ('GRAM' .in. uppercased) then
      ON_scheme%val = ON_scheme_values%Gram_Schmidt
    else if ('NO_O' .in. uppercased) then
      ON_scheme%val = ON_scheme_values%No_ON
    else
      call WarningMessage(2,'Invalid ORTH keyword')
      call Error(4)
      return
    end if
  end if
  if (Key('PERI')) write_orb_per_iter = .true.
  !---  Process NECI commands
  if (Key('NECI')) then
    if (DBG) write(u6,*) 'NECI is actived'
    DoNECI = .true.

    if (Key('DMPO')) then
      call WarningMessage(2,'NECI and DMPOnly are mutually exclusive.')
      call Error(4)
      return
    end if
    !-------------------------------------------------------------------
    if (Key('EMBD')) then
      DoEmbdNECI = .true.
#     ifndef _NECI_
      call WarningMessage(2,'EmbdNECI is given in input, so the embedded NECI should be used. Unfortunately MOLCAS was not '// &
                          'compiled with embedded NECI. Please use -DNECI=ON for compiling or use an external NECI.')
#     endif
    end if
    !-------------------------------------------------------------------
    if (Key('WRMA')) then
      WRMA = .true.
      if (DBG) write(u6,*) 'DMAT/PSMAT/PAMAT will be dumped.'
    end if
    !-------------------------------------------------------------------
    if (Key('GUGA')) then
      tGUGA_in = .true.
      if (DBG) write(u6,*) 'spin-free GUGA-NECI RDMs are actived'
      if (.not. Key('NECI')) then
        call WarningMessage(2,'GUGA requires NECI keyword!')
        call Error(4)
        return
      end if
    end if
    !--- This block is to process the DEFINEDET
    if (Key('DEFI')) then
      call setpos(luinput,'DEFI',line,irc)
      call mma_allocate(buffer,2000)
      if (irc /= _RC_ALL_IS_WELL_) then
        call Error(1)
        return
      end if
      ReadStatus = ' Failure reading Definedet.'
      read(luinput,'(A)',iostat=istatus) buffer
      if (istatus < 0) then
        call Error(2)
        return
      else if (istatus > 0) then
        call Error(3)
        return
      end if
      ReadStatus = ' O.K. reading Definedet.'
      call mma_allocate(definedet,len_trim(buffer))
      definedet(:) = trim(buffer)
      call mma_deallocate(buffer)
      write(u6,*) 'definedet read in proc_inp of size:',nactel
      write(u6,*) definedet
    end if
    if (Key('TOTA')) then
      call setpos(luinput,'TOTA',line,irc)
      if (irc /= _RC_ALL_IS_WELL_) then
        call Error(1)
        return
      end if
      ReadStatus = ' Failure reading totalwalkers.'
      read(luinput,*,iostat=istatus) totalwalkers
      if (istatus < 0) then
        call Error(2)
        return
      else if (istatus > 0) then
        call Error(3)
        return
      end if
      ReadStatus = ' O.K. reading totalwalkers.'
    else
      call WarningMessage(2,'TOTAlwalkers required for NECI.')
      call Error(4)
      return
    end if
    if (count([Key('RDML'),(Key('RDMS') .and. Key('CALC'))]) /= 1) then
      call WarningMessage(2,'RDMLinspace, and (RDMSamplingiters + CALCrdmonfly) are mutually exclusive, but one is required.')
      call Error(4)
      return
    else if (Key('RDML')) then
      call setpos(luinput,'RDML',line,irc)
      if (irc /= _RC_ALL_IS_WELL_) then
        call Error(1)
        return
      end if
      read(luinput,*,iostat=istatus) RDMsampling%start,RDMsampling%n_samples,RDMsampling%step
      if (istatus < 0) then
        call Error(2)
        return
      else if (istatus > 0) then
        call Error(3)
        return
      end if
    else if (Key('RDMS') .or. Key('CALC')) then
      if (.not. (Key('RDMS') .and. Key('CALC'))) then
        call WarningMessage(2,'RDMSamplingiters and CALCrdmonfly are both required.')
        call Error(4)
        return
      end if
      call setpos(luinput,'RDMS',line,irc)
      if (irc /= _RC_ALL_IS_WELL_) then
        call Error(1)
        return
      end if
      read(luinput,*,iostat=istatus) start
      if (istatus < 0) then
        call Error(2)
        return
      else if (istatus > 0) then
        call Error(3)
        return
      end if

      call setpos(luinput,'CALC',line,irc)
      if (irc /= _RC_ALL_IS_WELL_) then
        call Error(1)
        return
      end if
      read(luinput,*,iostat=istatus) length,step
      if (istatus < 0) then
        call Error(2)
        return
      else if (istatus > 0) then
        call Error(3)
        return
      end if
      RDMsampling = t_RDMsampling(start,length/step,step)
    end if
    if (Key('DIAG')) then
      call setpos(luinput,'DIAG',line,irc)
      if (irc /= _RC_ALL_IS_WELL_) then
        call Error(1)
        return
      end if
      read(luinput,*,iostat=istatus) diagshift
      if (istatus < 0) then
        call Error(2)
        return
      else if (istatus > 0) then
        call Error(3)
        return
      end if
    end if
    if (Key('TIME')) then
      call setpos(luinput,'TIME',line,irc)
      if (irc /= _RC_ALL_IS_WELL_) then
        call Error(1)
        return
      end if
      read(luinput,*,iostat=istatus) Time
      if (istatus < 0) then
        call Error(2)
        return
      else if (istatus > 0) then
        call Error(3)
        return
      end if
    end if
    if (Key('NMCY')) then
      call setpos(luinput,'NMCY',line,irc)
      if (irc /= _RC_ALL_IS_WELL_) then
        call Error(1)
        return
      end if
      read(luinput,*,iostat=istatus) nmCyc
      if (istatus < 0) then
        call Error(2)
        return
      else if (istatus > 0) then
        call Error(3)
        return
      end if
    end if
    if (Key('REAL')) then
      call setpos(luinput,'REAL',line,irc)
      if (irc /= _RC_ALL_IS_WELL_) then
        call Error(1)
        return
      end if
      read(luinput,*,iostat=istatus) realspawncutoff
      if (istatus < 0) then
        call Error(2)
        return
      else if (istatus > 0) then
        call Error(3)
        return
      end if
    end if
    if (Key('TRIA')) then
      call setpos(luinput,'TRIA',line,irc)
      if (irc /= _RC_ALL_IS_WELL_) then
        call Error(1)
        return
      end if
      read(luinput,*,iostat=istatus) trial_wavefunction
      if (istatus < 0) then
        call Error(2)
        return
      else if (istatus > 0) then
        call Error(3)
        return
      end if
    end if
    if (Key('POPS')) then
      call setpos(luinput,'POPS',line,irc)
      if (irc /= _RC_ALL_IS_WELL_) then
        call Error(1)
        return
      end if
      read(luinput,*,iostat=istatus) pops_trial
      if (istatus < 0) then
        call Error(2)
        return
      else if (istatus > 0) then
        call Error(3)
        return
      end if
    end if
    if (Key('SEMI')) then
      call setpos(luinput,'SEMI',line,irc)
      if (irc /= _RC_ALL_IS_WELL_) then
        call Error(1)
        return
      end if
      read(luinput,*,iostat=istatus) semi_stochastic
      if (istatus < 0) then
        call Error(2)
        return
      else if (istatus > 0) then
        call Error(3)
        return
      end if
    end if
    if (Key('MEMO')) then
      call setpos(luinput,'MEMO',line,irc)
      if (irc /= _RC_ALL_IS_WELL_) then
        call Error(1)
        return
      end if
      read(luinput,*,iostat=istatus) memoryfacspawn
      if (istatus < 0) then
        call Error(2)
        return
      else if (istatus > 0) then
        call Error(3)
        return
      end if
    end if
    !call fciqmc_option_check(iDoGas, nGSSH, iGSOCCX)
  end if

  if (Key('CCCI')) then
    if (DBG) write(u6,*) 'CC-CI is actived'
    Do_CC_CI = .true.

    if (Key('DMPO')) then
      call WarningMessage(2,'CC-CI and DMPOnly are mutually exclusive.')
      call Error(4)
      return
    end if
  end if

  ! ====================================================================
  !---  Process SYMM command
  if (Key('SYMM')) then
    if (DBG) write(u6,*) ' SYMMETRY command was given.'
    call SetPos(LUInput,'SYMM',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    Line = Get_Ln(LUInput)
    ReadStatus = ' Failure reading symmetry index after SYMM keyword.'
    read(Line,*,iostat=istatus) STSYM
    if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. reading symmetry index after SYMM keyword.'
    if (DBG) write(u6,*) ' State symmetry index ',STSYM
    call ChkIfKey()
  else
    ! If STSYM has not been set, normally it should be defaulted to 1.
    ! Exception: if this is a high-spin OS case, these often require STSYM /= 1:
    STSYM = 1
    if (ISPIN == NASHT+1) then
      do ISYM=1,NSYM
        NA = NASH(ISYM)
        if (mod(NA,2) /= 0) STSYM = MUL(STSYM,ISYM)
      end do
    end if
  end if
  call put_iscalar('STSYM',STSYM)
  if (DBG) write(u6,*) ' State symmetry STSYM=',STSYM

  ! ====================================================================

  !---  Process CIRE command
  if (Key('CIRE')) then
    if (DBG) write(u6,*) ' CIRESTART keyword was given.'
    ICIRST = 1
  end if

  !---  Process HOME command (root homing in SXCI part)
  if (Key('HOME')) then
    SXSEL = 'HOMING  '
    if (DBG) write(u6,*) ' HOME (Root Homing) keyword was given.'
    call SetPos(LUInput,'HOME',Line,iRc)
    call ChkIfKey()
  end if

  !---  Process SUPS command
  if (Key('SUPS')) then
    if (DBG) write(u6,*) ' SUPS (Supersymmetry) keyword was given.'
    call SetPos(LUInput,'SUPS',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    call mma_allocate(Temp1,mxOrb,Label='Temp1')
    ISUPSM = 1
    iOffset = 0
    do iSym=1,nSym
      ReadStatus = ' Failure reading data following SUPS keyword.'
      read(LUInput,*,iostat=istatus) nGrp
      if (istatus < 0) then
        call Error(2)
        return
      else if (istatus > 0) then
        call Error(3)
        return
      end if
      ReadStatus = ' O.K. after reading data following SUPS keyword.'
      do iGrp=1,nGrp
        call RdSups(LUInput,kOrb,Temp1)
        do iOrb=1,kOrb
          IXSYM(Temp1(iOrb)+iOffset) = iGrp
        end do
      end do
      iOffset = iOffset+nBas(iSym)
    end do
    call mma_deallocate(Temp1)
    ! (SVC) If both ALTER and SUPS keyword has been used, then change the IXSYM
    ! arrays according to the changed orbital ordering given in ALTER.
    do iAlter=1,NAlter
      iChng1 = IXSYM(iMAlter(iAlter,1))
      iChng2 = IXSYM(iMAlter(iAlter,2))
      IXSYM(iMAlter(iAlter,1)) = iChng2
      IXSYM(iMAlter(iAlter,2)) = iChng1
    end do
    call ChkIfKey()
  end if

  ! --- Process HEXS command
  if (Key('HEXS')) then
    if (DBG) write(u6,*) ' HEXS (Highly excited states) keyword was given. '
    call SetPos(LUInput,'HEXS',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    if ((I_ELIMINATE_GAS_MOLCAS /= 0) .and. (I_ELIMINATE_GAS_MOLCAS /= 2)) then
      call WarningMessage(2,'HEXS keyword defined more than once')
      call Error(1)
      return
    end if
    I_ELIMINATE_GAS_MOLCAS = I_ELIMINATE_GAS_MOLCAS+1
    ReadStatus = ' Failure reading data following HEXS keyword.'
    read(LUInput,*,iostat=istatus) N_ELIMINATED_GAS_MOLCAS
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading data following HEXS keyword.'
    ReadStatus = ' Failure reading data following HEXS keyword.'
    read(LUInput,*,iostat=istatus) (IELIMINATED_IN_GAS_MOLCAS(I),I=1,N_ELIMINATED_GAS_MOLCAS)
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading data following HEXS keyword.'
  end if

  ! --- Process DEXS command
  if (Key('DEXS')) then
    if (DBG) write(u6,*) ' DEXS (Doubly excited states) keyword was given. '
    call SetPos(LUInput,'DEXS',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    if (I_ELIMINATE_GAS_MOLCAS > 1) then
      call WarningMessage(2,'DEXS keyword defined more than once')
      call Error(1)
      return
    end if
    I_ELIMINATE_GAS_MOLCAS = I_ELIMINATE_GAS_MOLCAS+2
    ReadStatus = ' Failure reading data following DEXS keyword.'
    read(LUInput,*,iostat=istatus) N_2ELIMINATED_GAS_MOLCAS
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading data following DEXS keyword.'
    ReadStatus = ' Failure reading data following DEXS keyword.'
    read(LUInput,*,iostat=istatus) (I2ELIMINATED_IN_GAS_MOLCAS(I),I=1,N_2ELIMINATED_GAS_MOLCAS)
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading data following DEXS keyword.'
  end if

  !---  Process HROO command ---
  if (Key('HROO')) then
    if (DBG) write(u6,*) ' HROO (Hidden roots) keyword was given. '
    call SetPos(LUInput,'HROO',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading data following HROO keyword.'
    read(LUInput,*,iostat=istatus) hRoots
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading data following HROO keyword.'
  end if
  !---  Process NKEE command ---
  if (Key('NKEE')) then
    if (DBG) write(u6,*) ' NKEE (nr of kept vectors) keyword was given. '
    call SetPos(LUInput,'NKEE',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading data following NKEE keyword.'
    read(LUInput,*,iostat=istatus) n_keep
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading data following NKEE keyword.'
    if (n_keep < lRoots) then
      call WarningMessage(2,'nkeep must be at least equal to the number of roots')
      call Quit(_RC_INPUT_ERROR_)
    end if
  end if

  !---  Process CLEA command ---
  if (Key('CLEA')) then
    if (DBG) then
      write(u6,*) ' CLEAN (Orbital Cleaning) keyword.'
      write(u6,*) ' (Awkward input -- replace??).'
    end if
    call SetPos(LUInput,'CLEA',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    call mma_allocate(Temp1,mxOrb,Label='Temp1')
    call mma_allocate(Temp2,mxOrb,Label='Temp2')
    call mma_allocate(Temp3,mxOrb,Label='Temp3')
    nClean = sum(nBas(1:nSym)**2)
    call mma_allocate(Cleanmask,nClean,Label='CleanMask')
    CleanMask(:) = 0
    iOffset = 0
    do iSym=1,nSym
      mBas = nBas(iSym)
      ReadStatus = ' Failure reading data following CLEAN keyword.'
      read(LUInput,*,iostat=istatus) nGrp
      if (istatus < 0) then
        call Error(2)
        return
      else if (istatus > 0) then
        call Error(3)
        return
      end if
      ReadStatus = ' O.K. after reading data following CLEAN keyword.'
      do iGrp=1,nGrp
        call RdSups(LUInput,mOrb,Temp1)
        call RdSups(LUInput,nCof,Temp2)
        call RdSups(LUInput,mCof,Temp3)
        do i=1,mBas
          ii = (i-1)*mBas
          is_in_Group = 0
          do j=1,mOrb
            if (Temp1(j) == i) is_in_Group = 1
          end do
          if (is_in_Group == 1) then
            do k=1,nCof
              ij = Temp2(k)+ii+iOffset
              CleanMask(ij) = 1
            end do
          else
            do k=1,mCof
              ij = Temp3(k)+ii+iOffset
              CleanMask(ij) = 1
            end do
          end if
        end do
      end do
      iOffset = iOffset+mBas*mBas
    end do
    call mma_deallocate(Temp1)
    call mma_deallocate(Temp2)
    call mma_deallocate(Temp3)
    call ChkIfKey()
  end if

  !---  Process CHOL command (Cholesky Default Input, F.Aquilante Sept 04)
  if (Key('CHOL')) then
    if (DBG) write(u6,*) ' CHOLESKY keyword was given.'
    DoCholesky = .true.
    call SetPos(LUInput,'CHOL',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    call Cho_rasscf_rdInp(.true.,LUInput)
  end if

  !---  Process CHOI command (Cholesky Custom Input, F.Aquilante)
  ! Cholesky with user-defined settings.
  if (Key('CHOI')) then
    if (DBG) write(u6,*) ' CHOINPUT keyword was given.'
    DoCholesky = .true.
    call SetPos(LUInput,'CHOI',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    call Cho_rasscf_rdInp(.false.,LUInput)
  end if

  !---  Process ITER command
  if (Key('ITER')) then
    if (DBG) write(u6,*) ' ITERATIONS keyword was given.'
    call SetPos(LUInput,'ITER',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading data following ITER keyword.'
    read(LUInput,*,iostat=istatus) MAXIT,ITMAX
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. reading data following ITER keyword.'
    if (DBG) then
      write(u6,*) ' Max nr of RASSCF (macro) iterations MAXIT=',MAXIT
      write(u6,*) ' Max nr of orbital optimization iter ITMAX=',ITMAX
      if (Key('CION')) write(u6,*) ' (Irrelevant, since the CIONLY keyword was also given)'
    end if
    call ChkIfKey()
  end if

  !---  Process CRPR command
  if (Key('CRPR')) then
    if (DBG) write(u6,*) ' CRPR (Core Projector) keyword was used.'
    call SetPos(LUInput,'CRPR',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading core orbital after CRPR keyword.'
    read(LUInput,*,iostat=istatus) ITCORE,CORESHIFT
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. reading core orbital after CRPR keyword.'
    IfCRPR = .true.
    if (DBG) write(u6,*) ' Core orbital is number ITCORE'
    call ChkIfKey()
  end if

  !---  Process LEVS command
  if (Key('LEVS')) then
    if (DBG) write(u6,*) ' LEVS (Level Shift) keyword was used.'
    call SetPos(LUInput,'LEVS',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading level shift after LEVSHFT keyword.'
    read(LUInput,*,iostat=istatus) LVSHFT
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. reading level shift after LEVSHFT keyword.'
    if (DBG) write(u6,*) ' Level shift LVSHFT (Re*8!!) =',LVSHFT
    call ChkIfKey()
  end if

  !---  Process THRS command
  if (Key('THRS')) then
    if (DBG) write(u6,*) ' THRS (Thresholds) command was used.'
    call SetPos(LUInput,'THRS',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading thresholds after THRS keyword.'
    read(LUInput,*,iostat=istatus) THRE,THRTE,THRSX
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading thresholds after THRS keyword.'
    !tbp, may 2013: no altering of thre with Cholesky
    !tbp if (DoCholesky .and. (IPRLEV >= TERSE)) then
    !tbp   write(u6,*) '*** Detected Cholesky or RI/DF calculation'
    !tbp   write(u6,*) ' BUT user specified value of ThrE will be used. ThrE= ',THRE
    !tbp end if
    if (DBG) then
      write(u6,*) ' Threshold for energy change,    THRE =',THRE
      write(u6,*) ' Threshold for orbital rotation, THRTE=',THRTE
      write(u6,*) ' Threshold for max BLB element,  THRSX=',THRSX
    end if
    call ChkIfKey()
  end if

  !---  Process TIGH command
  if (Key('TIGH')) then
    if (DBG) write(u6,*) ' TIGHT (Tight CI convergence)  used.'
    KTIGHT = 1
    call SetPos(LUInput,'TIGH',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading data after TIGHT keyword.'
    read(LUInput,*,iostat=istatus) THREN,THFACT
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. reading data after TIGHT keyword.'
    if (DBG) then
      write(u6,*) ' CI energy threshold in 1st iter, THREN=',THREN
      write(u6,*) ' CI Threshold/Energy change =    THFACT=',THFACT
    end if
    call ChkIfKey()
  end if

  ! Use of quasi-newton ci/orbital coupling?
  ! Commands QUNE or NOQUNE:
  if (Key('NOQU')) then
    NQUNE = 0
    if (DBG) write(u6,*) ' NOQUNE keyword: QUNE is disabled.'
    call SetPos(LUInput,'NOQU',Line,iRc)
    call ChkIfKey()
  else if (Key('QUNE')) then
    if (DBG) write(u6,*) ' QUNE keyword: QUNE is enabled.'
    NQUNE = 1
    call SetPos(LUInput,'QUNE',Line,iRc)
    call ChkIfKey()
  else
    ! Default is to use QUNE, unless this is some kind of DFT:
    if (Key('FUNC')) then
      NQUNE = 0
      if (DBG) write(u6,*) ' DFT calculation: QUNE is disabled.'
    else
      NQUNE = 1
      if (DBG) write(u6,*) ' QUNE is enabled by default.'
    end if
  end if

  !---  Process CIMX command
  if (Key('CIMX')) then
    if (DBG) write(u6,*) ' Keyword CIMX was used.'
    call SetPos(LUInput,'CIMX',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading data after CIMX keyword.'
    read(LUInput,*,iostat=istatus) MAXJT
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading data after CIMX keyword.'
    if (DBG) write(u6,*) ' Max nr of CI iterations MAXJT=',MAXJT
    call ChkIfKey()
  end if

  !---  Process SDAV command
  if (Key('SDAV')) then
    if (DBG) write(u6,*) ' SDAV (Size of explicit Hamiltonian matrix)'
    call SetPos(LUInput,'SDAV',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading data after SDAV keyword.'
    read(LUInput,*,iostat=istatus) NSEL
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading data after SDAV keyword.'
    if (DBG) write(u6,*) ' Size is NSEL=',NSEL
    call ChkIfKey()
  end if

  !---  Process OFEM commands for Orbital-Free embedding
  if (Key('OFEM')) then
    if (DBG) write(u6,*) ' OFEM (Orbital-Free Embedding activated)'
    Do_OFemb = .true.
    call SetPos(LuInput,'OFEM',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading data after OFEM keyword.'
    read(LUInput,'(A)',iostat=istatus) OFE_KSDFT
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading data after OFEM keyword.'
    call UpCase(OFE_KSDFT)
    OFE_KSDFT = adjustl(OFE_KSDFT)
    write(u6,*)
    write(u6,*) '  --------------------------------------'
    write(u6,*) '   Orbital-Free Embedding Calculation'
    write(u6,*) '  --------------------------------------'
    if (OFE_KSDFT(1:4) == 'LDTF') then
      write(u6,*) '    T_nad potential   : Thomas-Fermi    '
    else
      write(u6,*) '    T_nad potential   : ',OFE_KSDFT(1:4)
    end if
    if (KEonly) then
      write(u6,*) '    Exc_nad potential :  None           '
    else
      write(u6,*) '    Exc_nad potential : ',OFE_KSDFT(6:10)
    end if
    write(u6,*) '  --------------------------------------'
    write(u6,*)
    DFTFOCK = 'SCF '
  end if
  if (Key('FTHA')) then
    call SetPos(LuInput,'FTHA',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading data after FTHA keyword.'
    read(LUInput,*,iostat=istatus) ThrFThaw
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading data after FTHA keyword.'
  end if
  if (Key('DFMD')) then
    call SetPos(LuInput,'DFMD',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading data after DFMD keyword.'
    read(LUInput,*,iostat=istatus) dFMD,Xsigma
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading data after DFMD keyword.'
    !write(u6,'(A,F6.3)') '  Fract. correl. potent. DFMD=',dFMD
    !write(u6,*) '  --------------------------------------'
    if (dFMD+Xsigma < Zero) then
      write(u6,*) ' *** Warning: arguments to DFMD must be nonnegative!'
      write(u6,*) ' ***          I will take their ABS !!! '
      dFMD = abs(dFMD)
      Xsigma = abs(Xsigma)
    end if
  end if

  !---  Process BKAP command for BK type of approximation
  !     (Giovanni Li Manni J.:GLMJ) Nov 2011
  if (Key('BKAP')) then
    DoBKAP = .true.
    call SetPos(LUInput,'BKAP',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading data after BKAP keyword.'
    read(LUInput,*,iostat=istatus) NGASBK
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    read(LUInput,*,iostat=istatus) (IOCCPSPC(IGAS,1),IGAS=1,NGASBK)
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    read(LUInput,*,iostat=istatus) (IOCCPSPC(IGAS,2),IGAS=1,NGASBK)
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. reading data after BKAP keyword.'
    if (DBG) then
      write(u6,*) ' BKAP: BK-type of approximation in action'
      write(u6,*) ' Min and Max for subspace with exact Hamiltonian '
      write(u6,*) ' =============================================== '
      write(u6,*) ' NGASBK :',NGASBK
      write(u6,*) '              Min. Occ.      Max. Occ.           '
      do IGAS=1,NGASBK
        write(u6,'(A,I2,10X,I3,9X,I3)') '   GAS',IGAS,IOCCPSPC(IGAS,1),IOCCPSPC(IGAS,2)
      end do
    end if
  end if

  !---  Process SPLI command for SplitCAS calculations
  !     (Giovanni Li Manni J.:GLMJ)
  if (Key('SPLI')) then
    if (DBG) write(u6,*) ' SPLI (Activation SplitCAS)'
    DoSplitCAS = .true.
    EnerSplit = .true.
    iDimBlockA = 0
    !*** The energy gap (GapSpli) is in eV ****
    GapSpli = Half*auToeV
    lrootSplit = 1
    thrSplit = 1.0e-6_wp
    MxIterSplit = 50
    call SetPos(LUInput,'SPLI',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
  end if

  !---  Process NUSP command for Numerical SplitCAS param. (GLMJ)
  if (Key('NUSP')) then
    if (DBG) write(u6,*) ' NUSP - Manual Setting of Numerical SplitCAS Param.'
    EnerSplit = .false.
    NumSplit = .true.
    call SetPos(LUInput,'NUSP',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading data after NUSP keyword.'
    read(LUInput,*,iostat=istatus) lrootSplit,iDimBlockA,MxIterSplit
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. reading data after NUSP keyword.'
    !read(LUInput,*,iostat=istatus) lrootSplit
    !if (istatus < 0) then
    !  call Error(2)
    !  return
    !else if (istatus > 0) then
    !  call Error(3)
    !  return
    !end if
    !read(LUInput,*,iostat=istatus) iDimBlockA
    !if (istatus < 0) then
    !  call Error(2)
    !  return
    !else if (istatus > 0) then
    !  call Error(3)
    !  return
    !end if
    !read(LUInput,*,iostat=istatus) MxIterSplit
    !if (istatus < 0) then
    !  call Error(2)
    !  return
    !else if (istatus > 0) then
    !  call Error(3)
    !  return
    !end if
    ReadStatus = ' Failure reading data after NUSP keyword.'
    read(LUInput,*,iostat=istatus) ThrSplit
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. reading data after NUSP keyword.'
    if (DBG) then
      write(u6,*) ' Root to be opt. in SplitCAS = ',lrootSplit
      write(u6,*) ' AA block size in SplitCAS = ',iDimBlockA
      write(u6,*) ' Max iteration in SplitCAS = ',MxIterSplit
      write(u6,*) ' Root to be opt. in SplitCAS = ',ThrSplit
    end if
  end if

  !---  Process ENSP command for Energetical SplitCAS param. (GLMJ)
  if (Key('ENSP')) then
    if (DBG) write(u6,*) ' ENSP - Manual Setting of Energetical SplitCAS Param.'
    EnerSplit = .true.
    NumSplit = .false.
    iDimBlockA = 0
    call SetPos(LUInput,'ENSP',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading data after ENSP keyword.'
    read(LUInput,*,iostat=istatus) lrootSplit,GapSpli,MxIterSplit
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. reading data after ENSP keyword.'
    !read(LUInput,*,iostat=istatus) lrootSplit
    !if (istatus < 0) then
    !  call Error(2)
    !  return
    !else if (istatus > 0) then
    !  call Error(3)
    !  return
    !end if
    !read(LUInput,*,iostat=istatus) GapSpli
    !if (istatus < 0) then
    !  call Error(2)
    !  return
    !else if (istatus > 0) then
    !  call Error(3)
    !  return
    !end if
    !read(LUInput,*,iostat=istatus) MxIterSplit
    !if (istatus < 0) then
    !  call Error(2)
    !  return
    !else if (istatus > 0) then
    !  call Error(3)
    !  return
    !end if
    ReadStatus = ' Failure reading data after ENSP keyword.'
    read(LUInput,*,iostat=istatus) ThrSplit
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. reading data after ENSP keyword.'
    if (DBG) then
      write(u6,*) ' Root to be opt. in SplitCAS = ',lrootSplit
      write(u6,*) ' Energy gap in SplitCAS = ',GapSpli
      write(u6,*) ' Max iteration in SplitCAS = ',MxIterSplit
      write(u6,*) ' Root to be opt. in SplitCAS = ',ThrSplit
    end if
  end if

  !---  Process PESP command for Percentage SplitCAS param. (GLMJ)
  if (Key('PESP')) then
    if (DBG) write(u6,*) ' PESP - Manual Setting of Percentage SplitCAS Param.'
    EnerSplit = .false.
    PerSplit = .true.
    NumSplit = .false.
    iDimBlockA = 0
    call SetPos(LUInput,'PESP',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading data after PESP keyword.'
    read(LUInput,*,iostat=istatus) lrootSplit,PercSpli,MxIterSplit
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    read(LUInput,*,iostat=istatus) ThrSplit
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. reading data after PESP keyword.'
    if (DBG) then
      write(u6,*) ' Root to be opt. in SplitCAS = ',lrootSplit
      write(u6,*) ' Percentage in SplitCAS = ',PercSpli
      write(u6,*) ' Max iteration in SplitCAS = ',MxIterSplit
      write(u6,*) ' Root to be opt. in SplitCAS = ',ThrSplit
    end if
  end if

  !---  Process FOSP command for First Order SplitCAS Approx. (GLMJ)
  if (Key('FOSP')) then
    if (DBG) write(u6,*) ' FOSP - First Order SplitCAS Approx.'
    FOrdSplit = .true.
  end if

  !---  Process OPTO keyword: Optimal Output for RASSCF/CASPT2
  !                           optimizations - GG Nov 2008.
  if (Key('OPTO')) then
    if (DBG) then
      write(u6,*) ' OPTO keyword was used.'
      write(u6,*) '(Optimal Output for RASSCF/CASPT2)'
    end if
    lOPTO = .true.
    call SetPos(LUInput,'OPTO',Line,iRc)
    call ChkIfKey()
  end if

  !---  Process SXDAmp command
  if (Key('SXDA')) then
    if (DBG) write(u6,*) ' SXDAMPING was requested.'
    call SetPos(LUInput,'SXDA',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading data after SXDAMP keyword.'
    read(LUInput,*,iostat=istatus) SXDamp
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading data after SXDAMP keyword.'
    if (DBG) write(u6,*) ' Parameter SXDamp=',SXDamp
    call ChkIfKey()
  end if

  !---  Process LOWM command
  if (Key('LOWM')) then
    if (DBG) then
      write(u6,*) ' LOWM keyword was used to force the CI routines'
      write(u6,*) ' to use Slater Determinants for low M and not M=S'
    end if
    LOWMS = 1
    call SetPos(LUInput,'LOWM',Line,iRc)
    call ChkIfKey()
  end if

  !---  Process LOWD keyword: Turn on Lowdin orthonormalization of CMOs
  if (Key('LOWD')) then
    if (DBG) then
      write(u6,*) ' LOWDIN orthonormalization was requested'
      write(u6,*) ' but from Jan 12 2010 that is default anyway.'
    end if
    Lowdin_ON = .true.
    call SetPos(LUInput,'LOWD',Line,iRc)
    call ChkIfKey()
  end if

  !---  Process PRWF command
  if (Key('PRWF')) then
    if (DBG) write(u6,*) ' The PRWF keyword was used.'
    call SetPos(LUInput,'PRWF',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading data after PRWF keyword.'
    read(LUInput,*,iostat=istatus) PRWTHR
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading data after PRWF keyword.'
    if (DBG) write(u6,*) ' Print CI coefficients larger than PRWTHR=',PRWTHR
    call ChkIfKey()
  end if

  !---  Process PRSD command
  if (Key('PRSD')) then
    if (DBG) write(u6,*) ' The PRSD keyword was used.'
    call SetPos(LUInput,'PRSD',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    if (DBG) write(u6,*) ' Print determinant expansions of CSFs'
    call ChkIfKey()
  end if

  !---  Process FCIDUMP command
  if (Key('FCID')) then
    ! activate the DMRG interface in RASSCF (dummy here since we stop after FCIDUMP)
    DOFCIDUMP = .true.
    call SetKey('DMRG',.true.)
    call SetPos(LUInput,'FCID',Line,iRc)
    call ChkIfKey()
  end if

# if ! defined (_ENABLE_BLOCK_DMRG_) && ! defined (_ENABLE_CHEMPS2_DMRG_)
  ! ====================================================================
  !          start of QCMaquis DMRG input section
  ! ====================================================================
# ifdef _DMRG_
  if (Key('DMRG') .or. doDMRG) then
    if (.not. doDMRG) then
      call SetPos(LUInput,'DMRG',Line,iRc)
      call ChkIfKey()
    end if
    !> DMRG flag
    doDMRG = .true.
    LRras2_dmrg(1:8) = 0
    !> LRras2 = Ras2 as the default
    LRras2_dmrg(1:nsym) = NRS2(1:nsym)
    !> initial guess setup
    guess_dmrg(1:7) = 'DEFAULT'
    call mma_allocate(initial_occ,nrs2t,nroots)
    initial_occ(:) = 0
  end if

  !---  Process RGIN command (QCMaquis Custom Input)
  if (Key('RGIN')) then
    if (DBG) write(u6,*) ' RGINPUT keyword was given.'
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    call SetPos(LUInput,'RGIN',Line,iRc)

    if (.not. Key('DMRG')) then
      if (.not. doDMRG) then
        call WarningMessage(2,'Error in input processing.')
        write(u6,*) ' PROC_INP: the keyword DMRG is not present but'
        write(u6,*) ' is required to enable the DMRG internal keyword'
        write(u6,*) ' section RGInput.'
        call Error(0)
        return
      end if
    end if

    nr_lines = 0
    call qcmaquis_rdinp(LuInput,1,nr_lines)
    call SetPos(LUInput,'RGIN',Line,iRc)
    call qcmaquis_rdinp(LuInput,2,nr_lines)

  end if

  !---  Process SOCC command (state occupation for initial guess in DMRG)
  if (Key('SOCC')) then
    if (DBG) write(u6,*) ' SOCC keyword was given.'
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    call SetPos(LUInput,'SOCC',Line,iRc)

    if (Key('DMRG') .or. doDMRG) then
      call socc_dmrg_rdinp(LuInput,initial_occ,nrs2t,nroots)
      guess_dmrg(1:7) = 'HF     '
    end if
  end if
# endif

  !-- Leon: Process NEVP(t2prep) keyword, prepare for 4-RDM calculation
  !    for (CD)-DMRG-NEVPT2
  if (Key('NEVP')) then
    if (DBG) write(u6,*) ' NEVP(t2prep) keyword was given.'
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
#   ifdef _DMRG_
    if ((.not. Key('DMRG')) .and. (.not. doDMRG)) then
      call WarningMessage(2,'Error in input processing.')
      write(u6,*) ' PROC_INP: the keyword DMRG is not present or'
      write(u6,*) ' DMRG not activated but is required to enable'
      write(u6,*) ' the NEVP keyword.'
      call Error(0)
      return
    end if

    DoNEVPT2Prep = .true.
    call SetPos(LUInput,'NEVP',Line,iRc)
    Line = Get_Ln(LUInput)
    call UpCase(Line)
    if (index(Line,'EVRD') /= 0) then
      call WarningMessage(2,'Warning,EvRDM keyword is deprecated.')
      write(u6,*) 'RDM evaluation is done in NEVPT2 module now'
    end if

#   else
    call WarningMessage(2,'Error in input processing.')
    write(u6,*) 'MOLCAS was compiled without QCMaquis support.'
    write(u6,*) 'Thus, no DMRG-NEVPT2 calculations are possible.'
    call Error(0)
    return
#   endif
  end if

# ifdef _DMRG_
  !> sanity checks
  !> a. DMRG requested but mandatory keywords not set at all
  if ((Key('DMRG') .or. doDMRG) .and. ((.not. Key('RGIN')) .and. (.not. as_solver_inp_proc))) then
    call WarningMessage(2,'Error in input processing.')
    write(u6,*) ' PROC_INP: the keyword RGINput is not present but'
    write(u6,*) ' is required for QCMaquis DMRG calculations in order'
    write(u6,*) ' to set compulsory DMRG internal parameters as'
    write(u6,*) ' for example:'
    write(u6,*) ' max_bond_dimension, conv_thresh and nsweeps.'
    write(u6,*) ' See the QCMaquis manual for further details.'
    call Error(0)
    return
  end if
  !> b. DMRG requested so check that ALL mandatory keywords have been set
  if ((Key('DMRG') .or. doDMRG) .and. (Key('RGIN') .or. as_solver_inp_proc)) then
    nr_lines = dmrg_input%nr_qcmaquis_input_lines
    call qcmaquis_rdinp(LuInput,3,nr_lines)
    if (nr_lines <= 0) then
      call Error(0)
      return
    end if
  end if
# endif
  ! ====================================================================
  !          end of QCMaquis DMRG input section
  ! ====================================================================
# endif

  !---  Process FARO command
  if (Key('FARO')) DoFaro = .true.

  !---  Process NOCA command
  if (DBG) write(u6,*) ' Check if NOCALC case.'
  if (Key('NOCA')) then
    if (DBG) write(u6,*) ' NOCALC keyword was used.'
    INOCALC = 1
    call SetPos(LUInput,'NOCA',Line,iRc)
    call ChkIfKey()
  end if

  !---  Process SAVE command
  if (DBG) write(u6,*) ' Check if SAVE_EXP case.'
  if (Key('SAVE')) then
    if (DBG) write(u6,*) ' SAVE_EXP keyword was used.'
    ISAVE_EXP = 1
    call SetPos(LUInput,'SAVE',Line,iRc)
    call ChkIfKey()
  end if

  !---  Process EXPA command
  if (DBG) write(u6,*) ' Check if EXPAND case.'
  if (Key('EXPA')) then
    if (DBG) write(u6,*) ' EXPAND keyword was used.'
    IEXPAND = 1
    call SetPos(LUInput,'EXPA',Line,iRc)
    call ChkIfKey()
  end if

  !---  Process DMRG command
# if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_)
  if (Key('DMRG')) then
    ! NN.14 FIXME: When DMRG option is disabled at compilation,
    !       this should give an error, but just ignored for the time.
    if (DBG) write(u6,*) ' DMRG (Use DMRG algorithm instead of FCI)'
    call SetPos(LUInput,'DMRG',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading data after DMRG keyword.'
    read(LUInput,*,iostat=istatus) MxDMRG
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading data after DMRG keyword.'
    if (DBG) write(u6,*) ' Nr. of states=',MxDMRG
    DoBlockDMRG = .true.
    call ChkIfKey()
  end if

  !---  Process 3RDM command
  if (Key('3RDM')) then
    if (DBG) write(u6,*) ' 3RDM (Compute 3RDM for DMRG-Cu4-CASPT2)'
    Do3RDM = .true.
#   ifdef _ENABLE_CHEMPS2_DMRG_
    iOrbTyp = 2
    IPT2 = 1
    write(u6,*) 'CHEMPS2> 3-RDM and F4-RDM require PseudoCanonical orbitals'
    write(u6,*) 'CHEMPS2> Automatically set: OUTOrbitals = CANOnical'
    if (Key('SUPS')) then
      write(u6,*) 'CHEMPS2> Bug using SYPSym and 3RDM, disable SUPSym'
      call Abend()
    end if
#   endif
    call SetPos(LUInput,'3RDM',Line,iRc)
    call ChkIfKey()
  end if

# ifdef _ENABLE_CHEMPS2_DMRG_
  !---  Process DAVT command
  if (Key('DAVT')) then
    call SetPos(LUInput,'DAVT',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading data after DAVT keyword.'
    read(LUInput,*,iostat=istatus) davidson_tol
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading data after DAVT keyword.'
    call ChkIfKey()
  end if

  !---  Process CHRE command
  if (Key('CHRE')) then
    if (DBG) write(u6,*) ' Restart in CheMPS2'
    chemps2_restart = .true.
    call SetPos(LUInput,'CHRE',Line,iRc)
    call ChkIfKey()
  end if

  !---  Process CHBL command
  if (Key('CHBL')) then
    call SetPos(LUInput,'CHBL',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading data after CHBL keyword.'
    read(LUInput,*,iostat=istatus) chemps2_blb
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading data after CHBL keyword.'
    call ChkIfKey()
  end if

  !---  Process MXSW command
  if (Key('MXSW')) then
    call SetPos(LUInput,'MXSW',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading data after MXSW keyword.'
    read(LUInput,*,iostat=istatus) max_sweep
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading data after MXSW keyword.'
    call ChkIfKey()
  end if

  !---  Process NOIS command
  if (Key('NOIS')) then
    call SetPos(LUInput,'NOIS',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading data after NOIS keyword.'
    read(LUInput,*,iostat=istatus) chemps2_noise
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading data after NOIS keyword.'
    call ChkIfKey()
  end if

  !---  Process DMRE command
  if (Key('DMRE')) then
    call SetPos(LUInput,'DMRE',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading data after DMRE keyword.'
    read(LUInput,*,iostat=istatus) chemps2_lrestart
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading data after DMRE keyword.'
    call ChkIfKey()
  end if

  !---  Process MXCA command
  if (Key('MXCA')) then
    call SetPos(LUInput,'MXCA',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading data after MXCA keyword.'
    read(LUInput,*,iostat=istatus) max_canonical
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading data after MXCA keyword.'
    call ChkIfKey()
  end if

# endif
# endif

  !---  Process HFOC command
  ! This keyword is to specify a user customized orbs occupancies guess.
  ! It is used by Block and CheMPS2... but it could be useful for other codes.
  ! Therefore it is now outside the ifdef Block or CheMPS2.
  if (Key('HFOC')) then
    write(u6,*) ' HFOC keyword was given.'
    call SetPos(LUInput,'HFOC',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading after HFOC keyword.'
    write(u6,*) 'NASHT, mxact = ',NASHT,mxact
    read(LUInput,*,iostat=istatus) (hfocc(i),i=1,NASHT)
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. reading after HFOC keyword.'
    write(u6,*) 'HFOCC read in proc_inp of size:',NASHT
    write(u6,*) (hfocc(i),i=1,NASHT)
  end if

# ifdef _ENABLE_DICE_SHCI_
  !---  Process DICE command
  if (Key('DICE')) then
    DoBlockDMRG = .true.
    write(u6,*) 'DICE> (semistochastic) heat bath configuration interaction (SHCI)'
    call SetPos(LUInput,'DICE',Line,iRc)
    call ChkIfKey()
  end if
  !---  Process STOC command
  if (Key('STOC')) then
    Dice_Stoc = .true.
    write(u6,*) 'DICE> Using semistochastic algorithm interaction (SHCI)'
    call SetPos(LUInput,'STOC',Line,iRc)
    call ChkIfKey()
  end if
  !---  Process DIOC command
  DICEOCC = ''
  if (Key('DIOC')) then
    call SetPos(LUInput,'DIOC',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading data after DIOC keyword.'
    read(LUInput,*,iostat=istatus) nref_dice
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    do iref_dice=1,nref_dice
      read(LUInput,'(A)',iostat=istatus) diceocc(iref_dice)
      if (istatus < 0) then
        call Error(2)
        return
      else if (istatus > 0) then
        call Error(3)
        return
      end if
      call molcas2dice(diceocc(iref_dice))
    end do
    ReadStatus = ' O.K. after reading data after DIOC keyword.'
    call ChkIfKey()
  end if
  !---  Process EPSI command
  if (Key('EPSI')) then
    if (DBG) write(u6,*) ' EPS (Thresholds) command was used.'
    call SetPos(LUInput,'EPS',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading thresholds after EPSI keyword.'
    read(LUInput,*,iostat=istatus) dice_eps1,dice_eps2
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading thresholds after EPSI keyword.'
    call ChkIfKey()
  end if
  !---  Process SAMP command
  if (Key('SAMP')) then
    call SetPos(LUInput,'SAMP',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading data after SAMP keyword.'
    read(LUInput,*,iostat=istatus) dice_sampleN
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading data after SAMP keyword.'
    call ChkIfKey()
  end if
  !---  Process DITE command
  if (Key('DITE')) then
    call SetPos(LUInput,'DITE',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) then
      call Error(1)
      return
    end if
    ReadStatus = ' Failure reading data after DITE keyword.'
    read(LUInput,*,iostat=istatus) dice_iter
    if (istatus < 0) then
      call Error(2)
      return
    else if (istatus > 0) then
      call Error(3)
      return
    end if
    ReadStatus = ' O.K. after reading data after DITE keyword.'
    call ChkIfKey()
  end if
  !---  Process DIRE command
  if (Key('DIRE')) then
    dice_restart = .true.
    call SetPos(LUInput,'DIRE',Line,iRc)
    call ChkIfKey()
  end if
# endif

  !---  All keywords have been processed ------------------------------*

  if ((.not. Key('INAC')) .and. (sum(nIsh(:nSym)) == 0)) &
    call WarningMessage(1,'The number of inactive orbitals is zero. Do you really want this?')

  if (IfCRPR) then
    ! Core shift using a fixed projection operator.
    NCRVEC = NBAS(1)
    call mma_allocate(CRVEC,NCRVEC,Label='CRVec')
    N = NBAS(1)
    NCRPROJ = (N*(N+1)*(N**2+N+2))/8
    call mma_allocate(CRPROJ,NCRPROJ,Label='CRPROJ')
  end if
  !*********************************************************************
  ! Generate artificial splitting or RAS into GAS for parallel blocking*
  !*********************************************************************
  if (.not. IDOGAS) then
    ! SVC: convert CAS/RAS to general GAS description here, then we only
    ! need to copy it for lucia later, which always uses GAS description.
    NGSSH(1,1:NSYM) = NRS1(1:NSYM)
    NGSSH(2,1:NSYM) = NRS2(1:NSYM)
    NGSSH(3,1:NSYM) = NRS3(1:NSYM)
    IGSOCCX(1,1) = max(2*sum(NRS1(1:NSYM))-NHOLE1,0)
    IGSOCCX(1,2) = 2*sum(NRS1(1:NSYM))
    IGSOCCX(2,1) = NACTEL-NELEC3
    IGSOCCX(2,2) = NACTEL
    IGSOCCX(3,1) = NACTEL
    IGSOCCX(3,2) = NACTEL
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Select default root for geometry optimization

  if ((NROOTS > 1) .and. (irlxroot == 0)) then

    ! Check if multi state SA-CASSCF

    nW = 0
    do iR=1,LROOTS
      if (WEIGHT(iR) /= Zero) nW = nW+1
    end do
    if (nW /= 1) then
      iRlxRoot = iroot(LROOTS)
    else
      do iR=1,LROOTS
        if (WEIGHT(iR) /= Zero) iRlxRoot = iroot(iR)
      end do
    end if
  end if
  if ((NROOTS == 1) .or. (LROOTS == 1)) iRlxRoot = iRoot(1)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Compute IZROT. IZROT is a matrix (lower triangular over the
  ! active space), which specifies which t,u rotations should be
  ! avoided, since the orbitals belong to the same RAS space.
  ! This is the only way the RAS concept is explicitly used in the
  ! SX section of the program.
  ITU = 0
  do ISYM=1,NSYM
    NAO = NASH(ISYM)

    if (DOBKAP) then
      !.Giovanni... BK stuff SplitCAS related. We want to treat RAS CI space as CAS.
      do NT=2,NAO
        IZROT(ITU+1:ITU+NT-1) = 1
        ITU = ITU+NT-1
      end do
    else

      do NT=2,NAO
        do NU=1,NT-1
          ITU = ITU+1
          IZROT(ITU) = 0
          !SVC: check if NU<NT are included in the same gas space
          NGSSH_LO = 0
          do IGAS=1,NGAS
            NGSSH_HI = NGSSH_LO+NGSSH(IGAS,ISYM)
            if ((NU > NGSSH_LO) .and. (NT <= NGSSH_HI)) IZROT(ITU) = 1
            NGSSH_LO = NGSSH_HI
          end do
        end do
      end do
    end if
  end do

  call Put_iArray('nIsh',nIsh,nSym)
  call Put_iArray('nAsh',nAsh,nSym)
  call Put_iScalar('Multiplicity',ISPIN)

  ! Initialize Cholesky information if requested
  if (DoCholesky) then
    call Cho_X_init(irc,ChFracMem)
    if (irc /= 0) then
      call Error(4)
      return
    end if
  end if

  ! ====================================================================

end if
!PAM Jump here in case of CASVB (IFVB == 2)
if (DBG) write(u6,*) ' After IFVB CONTINUE.'

!PAM July 2007 Check in case of CI restart:
if (DBG) write(u6,*) ' Check if CI-Restart.'
if (Key('CIRE')) then
  ! Test read:
  if (DBG) write(u6,*) ' Yes it is!'
  iJOB = -1
  call f_Inquire('JOBOLD',lExists)
  if (lExists) then
    if (DBG) write(u6,*) ' ''JOBOLD'' exists.'
    iJOB = 1
  else
    call f_Inquire(IPHNAME,lExists)
    if (lExists) then
      iJOB = 0
      if (DBG) then
        write(u6,*) ' No ''JOBOLD'', but JOBIPH exists.'
        write(u6,*) ' It is named ',IPHNAME
      end if
    end if
  end if
  if ((iJOB == 1) .or. (iJOB == 0)) then
    if (JOBOLD <= 0) JOBOLD = 20
    if (iJOB == 1) then
      call DaName(JOBOLD,'JOBOLD')
    else
      call DaName(JOBOLD,IPHNAME)
    end if
    IDisk = 0
    call IDaFile(JOBOLD,99,IScratch,10,IDisk)
    if ((JOBOLD > 0) .and. (JOBOLD /= JOBIPH)) then
      call DaClos(JOBOLD)
      JOBOLD = -1
    else if (JOBOLD > 0) then
      JOBOLD = -1
    end if
  else
    IScratch(1) = 0
  end if
  if (IScratch(1) == 0) then
    ! Test read failed. JOBOLD cannot be used.
    write(u6,*) ' Test read shows that there is no usable interface'
    write(u6,*) ' file, necessary for the requested CI restart.'
    write(u6,*) ' Most probable reason: the user has forgotten to'
    write(u6,*) ' provide this file. The program will continue,'
    write(u6,*) ' but there can be no CI restart.'
    ICIRST = 0
  end if
end if
!PAM July 2007 End of addition
! ======================================================================

! Initialize seward

if (DBG) write(u6,*) ' Initialize seward.'
nDiff = 0
if (DSCF .or. RF_On() .or. Langevin_On() .or. PCM_On() .or. Do_OFEmb .or. (KSDFT /= 'SCF')) &
  call IniSew(DSCF .or. Langevin_On() .or. PCM_On(),nDiff)
! ======================================================================
#ifdef _DMRG_
domcpdftDMRG = l_casdft .and. doDMRG
twordm_qcm = domcpdftDMRG .or. (.not. Key('CION'))
#endif

! Setup part for DMRG calculations
#ifdef _DMRG_
if (Key('DMRG') .or. doDMRG) then
  call getenvf('Project',ProjectName)
  call getenvf('WorkDir',WorkDir)
  ! Initialize the new interface

  ! in QCMaquis spins start with 0
  ! in QCMaquis irreps start with 0
  call qcmaquis_interface_init(nactel,sum(nrs2),ispin-1,stsym-1,nsym,nrs2,qcmaquis_param%conv_thresh,qcmaquis_param%M, &
                               qcmaquis_param%num_sweeps,trim(WorkDir)//'/'//trim(ProjectName),twordm_qcm,nroots,lroots,iroot, &
                               thre,weight, &
#ifdef _MOLCAS_MPP_
                               mpp_nprocs,mpp_procid, &
#endif
                               initial_occ)
  ! TODO: Support sweep_bond_dimension!
  ! This is an optional parameter to qcmaquis_interface_init not used here yet
  ! For now, it is set by the following qcmaquis_interface_set_param

  ! Read in the parameters
  ! Loop over the lines of qcmaquis input, ignore the comments
  do ii=1,size(dmrg_input%qcmaquis_input)/2
    ij = 2*ii
    call remove_comment(dmrg_input%qcmaquis_input(ij-1),'//')
    call remove_comment(dmrg_input%qcmaquis_input(ij),'//')

    call qcmaquis_interface_set_param(trim(dmrg_input%qcmaquis_input(ij-1)),trim(dmrg_input%qcmaquis_input(ij)))
  end do

  call qcmaquis_interface_stdout(trim(ProjectName)//'.QCMaquis.log')
end if
#endif

! Check the input data

if (DBG) write(u6,*) ' Call ChkInp.'
call ChkInp()
! ======================================================================

! In DMRG-CASSCF, skip GUGA and LUCIA settings
NCONF = 1
SkipGUGA = DoBlockDMRG
! ======================================================================

! Construct the Guga tables

if (.not. (DoNECI .or. Do_CC_CI .or. DumpOnly .or. SkipGUGA)) then
  ! right now skip most part of gugactl for GAS, but only call mknsm.
  if (.not. iDoGas) then
    ! DMRG calculation no need the GugaCtl subroutine
#   ifdef _DMRG_
    if (Key('DMRG') .or. doDMRG) then
      call mma_deallocate(initial_occ)
      SkipGUGA = .true.
    else
#   endif
      call Timing(Eterna_1,dum1,dum2,dum3)
      if (DBG) write(u6,*) ' Call GugaCtl'
      call GUGACtl(nSym,iSpin,nActEl,nHole1,nElec3,nRs1,nRs2,nRs3,SGS,CIS,EXS,STSYM,DoBlockDMRG)
      NCONF = CIS%NCSF(STSYM)

      call Timing(Eterna_2,dum1,dum2,dum3)
#   ifdef _DMRG_
    end if
#   endif
  else  ! if iDoGas
    call mknsm()
  end if
end if
! ======================================================================

if (.not. SkipGUGA) then
  ! Construct the determinant tables

  if (DBG) write(u6,*) ' Construct the determinant tables.'
  MS2 = iSpin-1

  ! Set variables needed in Lucia_Ini

  ngssh_Molcas(:,:) = ngssh(:,:)
  igsoccx_Molcas(:,:) = igsoccx(:,:)
  potnuc_Molcas = potnuc
  thre_Molcas = thre
  nsym_Molcas = nsym
  nactel_Molcas = nactel
  ms2_Molcas = ms2
  ispin_Molcas = ispin
  lsym_Molcas = stsym
  itmax_Molcas = itmax
  nroots_Molcas = max(nroots,lRoots)
  ipt2_Molcas = ipt2
  iprci_molcas = iprloc(3)
  ngas_molcas = ngas
  INOCALC_MOLCAS = INOCALC
  ISAVE_EXP_MOLCAS = ISAVE_EXP
  IEXPAND_MOLCAS = IEXPAND

  ! And call Lucia_Ini to initialize LUCIA

  ! Combinations don't work for CASVB (at least yet)!
  if (ifvb /= 0) iSpeed(1) = 0

  if (.not. (Key('DMRG') .or. DoNECI .or. Do_CC_CI .or. DumpOnly)) then
    ! switch on/off determinants
#   ifdef _DMRG_
    if (.not. doDMRG) then
#   endif
      ! Initialize LUCIA and determinant control
      call StatusLine('RASSCF: ','Initializing Lucia...')
      call Lucia_Util('Ini')
      ! to get number of CSFs for GAS
      ! and number of determinants to store
      nconf = sum(ncsasm(1:mxsym))
      nDet = sum(ndtasm(1:mxsym))
#   ifdef _DMRG_
    end if
#   endif
  end if

  ISCF = 0
  if ((ISPIN == NAC+1) .and. (NACTEL == NAC)) ISCF = 1
  if ((ISPIN == 1) .and. (NACTEL == 2*NAC)) ISCF = 1
  if (ISCF == 1) then
    NCONF = 1
    MAXJT = 1
  end if

  ! If the CI-root selectioning option has been specified translate
  ! the reference configuration numbers from the split graph GUGA
  ! to the symmetric group numbering

  ! ====================================================================
  if (ICICH == 1) then
    call mma_allocate(UG2SG_X,NCONF,Label='UG2SG_X')
    call UG2SG(NROOTS,NCONF,NAC,NACTEL,STSYM,IPR,CONF,CFTP,UG2SG_X,ICI,JCJ,CCI,MXROOT)
    call mma_deallocate(UG2SG_X)
  end if
  ! ====================================================================

  ! faroald initializations
  if (DOFARO) then
    if (NSYM > 1) then
      write(u6,'(1X,A)') 'FARO keyword was used, but NSYM > 1,'
      write(u6,'(1X,A)') 'switching to LUCIA as the CI backend.'
      DOFARO = .false.
    else
      write(u6,'(1X,A)') '**EXPERIMENTAL**'
      write(u6,'(1X,A)') 'CI backend is FAROALD instead of LUCIA.'
      write(u6,'(1X,A)') '**EXPERIMENTAL**'
      call FAROALD_INIT(NACTEL,NASH(1),ISPIN)
      call CITRANS_INIT(NACTEL,NASH(1),ISPIN)
    end if
  end if

end if

!---  Normal exit -----------------------------------------------------*
if (DBG) write(u6,*) ' Normal exit from PROC_INP.'

return

contains

subroutine Error(code)

  integer(kind=iwp), intent(in) :: code

  select case (code)
    case (0)
      iRc = _RC_INPUT_ERROR_
    case (1)
      if (IPRLEV >= TERSE) then
        call WarningMessage(2,'Error in input preprocessing.')
        write(u6,*) ' PROC_INP: A keyword was found during prescanning'
        write(u6,*) ' the input file, but when later trying to locate'
        write(u6,*) ' this input, it could not be found. Something has'
        write(u6,*) ' happened to the input file, or else there is some'
        write(u6,*) ' strange program error.'
      end if
    case (2)
      call WarningMessage(2,'End of input file during preprocessing.')
      call WarningMessage(2,ReadStatus)
      if (IPRLEV >= TERSE) write(u6,*) ' Error exit from PROC_INP.'
    case (3)
      call WarningMessage(2,'Read error during input preprocessing.')
      call WarningMessage(2,ReadStatus)
      if (IPRLEV >= TERSE) write(u6,*) ' Error exit from PROC_INP.'
    case (4)
      call WarningMessage(2,'Error during input preprocessing.')
      call WarningMessage(2,ReadStatus)
      if (IPRLEV >= TERSE) then
        write(u6,*) ' Error exit from PROC_INP.'
        write(u6,*) ' Check previous messages in the output'
        write(u6,*) ' to find the reason.'
      end if
      iRc = _RC_INPUT_ERROR_
  end select
  if (DBG) write(u6,*) ' Abnormal exit from PROC_INP.'

end subroutine Error

end subroutine proc_inp
