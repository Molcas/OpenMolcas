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

use fortran_strings, only: to_upper, operator(.in.)
use csfbas, only: CONF
use lucia_data, only: CFTP
use Fock_util_global, only: DoCholesky
use Cholesky, only: ChFracMem
use write_orbital_files, only: OrbFiles, write_orb_per_iter
use fcidump, only: DumpOnly
use fcidump_reorder, only: ReOrInp, ReOrFlag
use fciqmc, only: DoEmbdNECI, DoNECI, tGUGA_in, tPrepStochCASPT2, tNonDiagStochPT2
use fciqmc_read_RDM, only: MCM7, WRMA
use CC_CI_mod, only: Do_CC_CI
use spin_correlation, only: orb_range_p, orb_range_q, same_orbs, tRootGrad
use orthonormalization, only: ON_scheme, ON_scheme_values
use fciqmc_make_inp, only: trial_wavefunction, pops_trial, t_RDMsampling, RDMsampling, totalwalkers, Time, nmCyc, memoryfacspawn, &
                           realspawncutoff, diagshift, definedet, semi_stochastic
use casvb_global, only: ifvb
use KSDFT_Info, only: CoefR, CoefX
use OFembed, only: Do_OFemb, KEonly, OFE_KSDFT, ThrFThaw, Xsigma, dFMD
use CMS, only: iCMSOpt, CMSGiveOpt, CMSGuessFile
use UnixInfo, only: SuperName
use Lucia_Interface, only: Lucia_Util
use gugx, only: SGS, CIS, EXS
use gas_data, only: iDoGAS, NGAS, NGSSH, IGSOCCX
use Symmetry_info, only: Mul
use SplitCas_Data, only: DoSPlitCas, MxIterSplit, ThrSplit, lRootSplit, NumSplit, EnerSplit, PerSplit, PerCSpli, fOrdSplit, &
                         iDimBlockA, GapSpli
use PrintLevel, only: DEBUG, VERBOSE, TERSE
use output_ras, only: LF, IPRGLB, IPRLOC
use general_data, only: MAXALTER, NALTER, JOBIPH, NSYM, INVEC, STARTORBFILE, NBAS, LUSTARTORB, JOBOLD, NTOT, NTOT1, NTOT2, NDELT, &
                        NFROT, NTOTSP, NRS1T, NRS2T, NRS3T, NACTEL, NHOLE1, NELEC3, ISPIN, STSYM, NSEL, SXDAMP, LOWDIN_ON, NISH, &
                        NCRVEC, NRS1, NRS2, NRS3, NCONF, MALTER, NASH, NDEL, NFRO, NORB, NSSH, CRVec, CleanMask, CRPROJ
use spinfo, only: NSYM_MOLCAS, NACTEL_MOLCAS, MS2_MOLCAS, ISPIN_MOLCAS, LSYM_MOLCAS, NROOTS_MOLCAS, NGAS_MOLCAS, THRE_MOLCAS, &
                  ITMAX_MOLCAS, INOCALC_MOLCAS, ISAVE_EXP_MOLCAS, IEXPAND_MOLCAS, IPT2_MOLCAS, I_ELIMINATE_GAS_MOLCAS, &
                  N_ELIMINATED_GAS_MOLCAS, N_2ELIMINATED_GAS_MOLCAS, IPRCI_MOLCAS, POTNUC_MOLCAS, I2ELIMINATED_IN_GAS_MOLCAS, &
                  IELIMINATED_IN_GAS_MOLCAS, IGSOCCX_MOLCAS, ISPEED, NGSSH_MOLCAS, DOBKAP, NGASBK, IOCCPSPC, MS2, NDET, NCSASM, &
                  NDTASM
use DWSol, only: DWSol_DWRO
use Molcas, only: LenIn, MxAct, MxGAS, MxOrb, MxRoot, MxSym
use RASDim, only: MxRef, MxTit
use input_ras   ! It should be without the only option!
use rasscf_global, only: KSDFT, IROOT, IRLXROOT, ICI, CCI, HFOCC, CMSStartMat, CMSThreshold, CoreShift, DFTFOCK, DoBLOCKDMRG, &
                         ExFac, hRoots, iAlphaBeta, ICICH, iCIonly, iCIRFROOT, iCMSITERMAX, iCMSITERMin, iCMSP, iExpand, JCJ, &
                         iFORDE, iOrbOnly, iOrbTyp, iOrdEM, iOverWr, iPCMRoot, iPhName, iPR, iPT2, iRotPsi, iSave_Exp, iSCF, &
                         iSPDEN, iSupSM, ITCORE, ITMAX, IXMSP, kivo, kTight, l_CASDFT, LowMS, LvShft, MaxIt, MaxJt, MaxOrbOut, &
                         n_keep, NAC, NACPAR, NACPR2, NFR, NIN, NO2M, NonEq, NQUNE, NORBT, NROOTS, NSEC, NTOT3, NTOT4, OutFmt1, &
                         OutFmt2, PotNuc, PreThr, ProThr, PreThr, Purify, RFPert, S, SXSEL, ThFact, ThrE, ThrEn, ThrSX, ThrTE, &
                         HFOcc, Title, Weight, DoFaro, DoFCIDump, iCIRST, IfCRPR, LROOTS, PrwThr, InOCalc, ixSym, iZRot
#ifdef _DMRG_
#endif
#ifdef _ENABLE_DICE_SHCI_
use rasscf_global, only: diceOcc, dice_eps1, dice_eps2, dice_iter, dice_restart, dice_sampleN, dice_stoc, nRef_dice
#endif
#ifdef _ENABLE_CHEMPS2_DMRG_
use rasscf_global, only: ChemPS2_Restart, ChemPS2_lRestart, Davidson_Tol, ChemPS2_BLB, Max_Sweep, ChemPS2_Noise, Max_Canonical, &
                         MxDMRG, Do3RDM
#endif
#ifdef _DMRG_
use qcmaquis_interface_cfg
use qcmaquis_interface, only: qcmaquis_interface_init, remove_comment, qcmaquis_interface_set_param, qcmaquis_interface_stdout
use active_space_solver_cfg, only: as_solver_inp_proc
use rasscf_global, only: MPSCompressM, DoNEVPT2Prep, Twordm_qcm, DoMCPDFTDMRG, DoDMRG
#ifdef _MOLCAS_MPP_
use Para_Info, only: mpp_procid, mpp_nprocs
#endif
#endif
#ifdef _HDF5_
use mh5, only: mh5_is_hdf5, mh5_open_file_r, mh5_exists_attr, mh5_exists_dset, mh5_fetch_attr, mh5_fetch_dset, mh5_close_file
#endif
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
logical lOPTO
logical DSCF
integer iRC
character(len=180) Line
character(len=8) NewJobIphName
logical lExists, RunFile_Exists, RlxRCheck
logical RF_On
logical Langevin_On
logical PCM_On
integer iMAlter(MaxAlter,2)
integer IPRGLB_IN, IPRLOC_IN(7)
logical DBG, exist
integer IScratch(10)
character(len=8) InfoLbl
integer NBAS_L(8), NORB_L(8)
integer NFRO_L(8), NISH_L(8), NRS1_L(8), NRS2_L(8)
integer NRS3_L(8), NSSH_L(8), NDEL_L(8)
integer IADR19(15)
character(len=180) Get_LN
external Get_LN
real*8 Get_ExFac
external Get_ExFac
character(len=72) ReadStatus
character(len=72) JobTit(mxTit)
character(len=256) myTitle
character(len=8) MaxLab
logical, external :: Is_First_Iter
real*8 Dummy(1)
character(len=(LenIn+8)*mxOrb) lJobH1
character(len=2*72) lJobH2
integer :: start, step, length
character(len=50) :: ON_scheme_inp, uppercased
character(len=:), allocatable :: buffer
intrinsic INDEX, NINT, DBLE, SQRT
integer, allocatable :: Temp1(:), Temp2(:), Temp3(:), type(:), Stab(:), UG2SG_X(:)
real*8, allocatable :: ENC(:), RF(:)
real*8 dSum, dum1, dum2, dum3, Eterna_2, POTNUCDUMMY, PRO, SUHF, TEffNChrg, TotChrg, Eterna_1
integer, external :: IsFreeUnit, nToken
integer i, i1, i2, iad19, iChng1, iChng2, iDisk, iEnd, iErr, iGAS, iGrp, ii, ij, iJOB, inporb_version, iod_save, iOffSet, iOrb, &
        iOrbData, iPrLev, iR, iRC1, iRef, iReturn, is_in_group, iStart, iSum, iSym, itu, j, jpcmroot, k, korb, kref, mBas, mCof, &
        mConf, mm, mOrb, N, NA, NAO, NASHT, NCHRG, nClean, nCof, nDiff, nGrp, NGSSH_HI, NGSSH_LO, NISHT, nItems, nNUc, nOrbRoot, &
        nOrbs, nSym_l, nT, nU, nW, iAll, iAlter, NISHT_old, NCRPROJ
#ifdef _ENABLE_DICE_SHCI_
integer iref_dice
#endif
#ifdef _DMRG_
character(len=256) WorkDir
character(len=72) ProjectName
integer :: LRras2_dmrg(8)
integer, allocatable :: initial_occ(:,:)
character(len=20) :: guess_dmrg
integer nr_lines
#endif
#ifdef _HDF5_
character(len=1), allocatable :: typestring(:)
integer mh5id, lRoots_l
#endif
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
davidson_tol = 1.0d-7
chemps2_blb = 0.5d-2
max_sweep = 8
chemps2_noise = 0.05
max_canonical = max_sweep*5
#endif
! Init HFOCC array containing user defined occupancies for the active orbitals.
! This array is used by DMRG codes (Block as well as CheMPS2).
! Terefore I took it out of any ifdef preprocessing flag.

do i=1,MxAct
  hfocc(i) = 0
end do

#ifdef _ENABLE_DICE_SHCI_
dice_stoc = .false.
nref_dice = 1
dice_eps1 = 1.0d-4
dice_eps2 = 1.0d-5
dice_sampleN = 200
dice_iter = 20
dice_restart = .false.
#endif

!    SplitCAS related variables declaration  (GLMJ)
DoSplitCAS = .false.
NumSplit = .false.
EnerSplit = .false.
PerSplit = .false.
FOrdSplit = .false.
!    BK type of approximation (GLMJ)
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
if (KeyCORE) then
  call SetPos(LUInput,'CORE',Line,iRc)
  call ChkIfKey()
end if
if (KeyINPO) then
  call SetPos(LUInput,'INPO',Line,iRc)
  call ChkIfKey()
end if
if (KeyJOBI) then
  call SetPos(LUInput,'JOBI',Line,iRc)
  call ChkIfKey()
end if
if (KeyLUMO) then
  call SetPos(LUInput,'LUMO',Line,iRc)
  call ChkIfKey()
end if
if (KeyCIRE) then
  call SetPos(LUInput,'CIRE',Line,iRc)
  call ChkIfKey()
end if
! They had to be tested here, since the key flags may be modified later.

! ======================================================================
! The outdated INPORB keyword now means just the same as LUMORB
! Will be deleted soon.
KeyLUMO = KEYLUMO .or. KeyINPO
! ======================================================================

! How was the program called?
!PAM 2009  For some particular types of calculations, the input is to be
! disregarded or overridden, as follows (Roland, esp. numerical diff):
if (KeyEXPE) then
  call SetPos(LUInput,'EXPE',Line,iRc)
  call ChkIfKey()
else
  IfVB = 0
  if (SuperName(1:6) == 'rasscf') then
    ! For geometry optimizations use the old CI coefficients.
    if (.not. Is_First_Iter()) then
      KeyCIRE = .true.
      KeyFILE = .false.
    end if
  else if (SuperName(1:5) == 'casvb') then
    IfVB = 2
  else if (SuperName(1:6) == 'loprop') then
    KeyCIRE = .true.
    KeyFILE = .false.
  else if (SuperName(1:11) == 'last_energy') then
    KeyCIRE = .true.
    KeyFILE = .false.
  else if (SuperName(1:18) == 'numerical_gradient') then
    KeyCIRE = .true.
    KeyFILE = .false.
  end if
end if
!PAM2009 Also, CIRESTART would normally also imply that orbitals are
!        to be taken from JOBIPH or JOBOLD:
if (.not. KeyEXPE) then
  if (KeyCIRE) then
    KeyLUMO = .false.
    KeyINPO = .false.
    KeyJOBI = .true.
  end if
end if

! Check PRINT command
if (KeyPRIN) then
  call SetPos(LUInput,'PRIN',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading after PRINT keyword.'
  read(LUInput,*,end=9910,err=9920) LF,IPRGLB_IN,(IPRLOC_IN(I),I=1,7)
  ReadStatus = ' O.K. reading info of PRINT keyword.'
  call SetPrLev(IPRGLB_IN,IPRLOC_IN)
  if (IPRLOC(1) >= DEBUG) then
    write(6,*) ' PROC_INP: Print levels have been set:'
    write(6,*) '    Global print level IPRGLB=',IPRGLB
    write(6,*) '    Local print levels by section:'
    write(6,*) '             Input section, IPRLOC(1)=',IPRLOC(1)
    write(6,*) '    Transformation section, IPRLOC(2)=',IPRLOC(2)
    write(6,*) '                CI section, IPRLOC(3)=',IPRLOC(3)
    write(6,*) '          Super-CI section, IPRLOC(4)=',IPRLOC(4)
    write(6,*) '            Output section, IPRLOC(6)=',IPRLOC(6)
    write(6,*) '          Property section, IPRLOC(7)=',IPRLOC(7)
  end if
  call ChkIfKey()
end if

! Local print level in this routine:
IPRLEV = IPRLOC(1)
! Short, for use with IF-statements at debugging print level:
DBG = DBG .or. (IPRLEV >= DEBUG)

if (DBG) write(6,*) ' Trace of input processing in PROC_INP:'

! ======================================================================
! Check IfVB flag
! if (IFVB == 2) -- Then bypass a lot of processing:
if (DBG) write(6,*) ' Valence-Bond flag IfVB=',IfVB
if (ifvb == 2) then
  if (JOBIPH <= 0) then
    JOBIPH = IsFreeUnit(15)
    call daname(JOBIPH,'JOBIPH')
  end if
  ! Special input routine:
  if (DBG) write(6,*) ' IfVB=2, so Call Readin_VB'
  call Readin_vb()
  if (DBG) write(6,*) ' Back from Readin_VB'
  if (DBG) write(6,*) ' Bypass usual input processing, GoTo 100!'
end if
if (ifvb == 2) goto 100

! ======================================================================
if (DBG) write(6,*) ' Check if VB keyword was given.'
!---  Process vb   command
if (KeyVB) then
  if (DBG) write(6,*) ' Yes it was!'
  call SetPos(LUInput,'VB  ',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  if (DBG) write(6,*) ' so Call CVBInp_CVB.'
  call cvbinp_rvb(1,LUInput)
  if (DoCholesky) then
    call WarningMessage(2,'CASVB cannot do Cholesky or RI/DF.')
    call Quit(_RC_INPUT_ERROR_)
  end if
  if (DBG) write(6,*) ' Set IfVB=1.'
  ifvb = 1
  if (DBG) write(6,*) ' Asked for CASVB calc.'
end if

! ======================================================================
!---  process TITLE    command
if (KeyTITL) then
  call SetPos(LUInput,'TITL',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  Title(1) = trim(Get_Ln(LUInput))
  if (DBG) then
    write(6,*) ' PROC_INP: Title line:'
    write(6,*) Title(1)
  end if
  call ChkIfKey()
end if

! ======================================================================
! Check if there is any runfile:
call F_Inquire('RUNFILE',RunFile_Exists)
if (DBG) write(6,*) ' Inquire about RUNFILE.'
if (RunFile_Exists) then
  if (DBG) write(6,*) ' Yes, there is one.'
  NSYM = 0
  call qpg_iScalar('nSym',lExists)
  if (lExists) then
    call Get_iScalar('nSym',nSym)
    call Get_iArray('nBas',nBas,nSym)
    if (DBG) then
      write(6,*) ' The following information exists on runfile:'
      write(6,*) ' Nr of symmetries, NSYM:',NSYM
      write(6,*) ' Nr of basis functions/symmetry:'
      write(6,'(1x,8I5)') (NBAS(I),I=1,NSYM)
      call XFlush(6)
    end if
  else
    call WarningMessage(2,'No symmetry info on runfile.')
    write(6,*) ' There seems to be no information about symmetry'
    write(6,*) ' on the runfile! This is an unexpected error.'
    call Quit(_RC_IO_ERROR_READ_)
  end if
else
  call WarningMessage(2,'Cannot find runfile.')
  write(6,*) ' PROC_INP: Cannot find RUNFILE. This is an unexpected error.'
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
if (DBG) write(6,*) ' Where to read orbitals? '
StartOrbFile = 'INPORB'
if (KeyFILE) then
  if (DBG) write(6,*) ' Reading file name for start orbitals.'
  call SetPos(LUInput,'FILE',Line,iRc)
  Line = Get_Ln(LUInput)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  call ChkIfKey()
  if (DBG) then
    write(6,*) ' Calling fileorb with filename='
    write(6,*) Line
  end if
  call fileorb(Line,StartOrbFile)
# ifdef _HDF5_
  if (mh5_is_hdf5(StartOrbFile)) then
    KeyLUMO = .false.
    KeyTYPE = .false.
    KeyH5OR = .true.
    if (KeyCIRE) then
      write(6,*) ' CIRE keyword given, but FILE is an HDF5 format.'
      write(6,*) ' I will switch off CIRE and use H5CI instead!'
      KeyH5CI = .true.
      KeyCIRE = .false.
      KeyJOBI = .false.
    end if
  else
    if (KeyCIRE) then
      write(6,*) ' CIRE keyword given, FILE will be ignored!'
    else
      KeyLUMO = .true.
    end if
  end if
# else
  if (KeyCIRE) then
    write(6,*) ' CIRE keyword given, FILE will be ignored!'
  else
    KeyLUMO = .true.
  end if
# endif
end if
if (DBG) then
  write(6,*) ' StartOrbFile='
  write(6,*) StartOrbFile
  call xflush(6)
end if

! ======================================================================
! The JOBIPH file, following decisions from the Torre Normanna labour camp:
! Default, the file name is 'JOBIPH'.
! However, if keyword IPHNAME was used, then the name was given in input.
! Also, if instead the keyword NEWIPH was given, then a new name will be
! chosen as the first not-already-used name in the sequence
! 'JOBIPH', 'JOBIPH01', 'JOBIPH02', etc.
if (DBG) write(6,*) ' Present name of JOBIPH file is ',IPHNAME
!---  Process NEWIPH command (P A Malmqvist Sep 06)
if (KeyNEWI) then
  if (DBG) write(6,*) ' A fresh JOBIPH file will be used.'
  IPHNAME = 'ToBeFoun'
end if
!---  process IPHNAME command (P A Malmqvist Sep 06)
if (KeyIPHN) then
  if (DBG) write(6,*) ' Reading file name for JOBIPH file.'
  call SetPos(LUInput,'IPHN',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading IPHNAME string.'
  read(LUInput,*,end=9910,err=9920) IPHNAME
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
      if (.not. lExists) goto 12
    end do
    write(LF,*)
    write(LF,*) '******************************************'
    write(LF,*) ' Sorry, all possible JOBIPH names in this '
    write(LF,*) ' directory already in use. Program stops. '
    write(LF,*) '******************************************'
    call Abend()
12  continue
    IPHNAME = NewJobIphName
  end if
end if
if (DBG) then
  write(6,*) ' Name of JOBIPH file is'
  write(6,*) IPHNAME
end if
! Finally, we know the name. Open jobiph file. If another file is already
! opened with this identifier close it.
if (JOBIPH > 0) then
  call DaClos(JOBIPH)
  JOBIPH = -1
end if
JOBIPH = IsFreeUnit(15)
call DANAME(JOBIPH,IPHNAME)

! ======================================================================
! If orbital files should be produced, several keywords are relevant:
!---  Process OUTPRINT command
if (KeyOUTP) then
  if (DBG) write(6,*) ' OUTPRINT command was given:'
  call SetPos(LUInput,'OUTP',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading OUTPRINT line.'
  read(LUInput,*,end=9910,err=9920) Line
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
    write(LF,*) ' Input line is:'
    write(LF,*) Line
    write(LF,*) ' Did not understand ''OUTP'' input. Ignored.'
  end if
end if
!---  Process PROR (Print levels for orbitals) command
if (KeyPROR) then
  if (DBG) write(6,*) ' PRORB command was given:'
  call SetPos(LUInput,'PROR',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading PRORB input.'
  read(LUInput,*,end=9910,err=9920) PRETHR,PRO
  ReadStatus = ' O.K. after reading PRORB input.'
  PROTHR = max(0.0d0,PRO)
end if
!---  Process ORBL keyword: Orbital listing
if (KeyORBL) then
  if (DBG) write(6,*) ' ORBL (Orbital listing):'
  call SetPos(LUInput,'ORBL',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading ORBL input.'
  read(LUInput,*,end=9910,err=9920) Line
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
    write(LF,*) ' Input line is:'
    write(LF,*) Line
    write(LF,*) ' Did not understand ''OUTL'' input. Ignored.'
  end if
end if
!---  Process ORBA keyword: Orbital Appearance
if (KeyORBA) then
  if (DBG) write(6,*) ' ORBA (Orbital appearance):'
  call SetPos(LUInput,'ORBA',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading ORBA input.'
  read(LUInput,*,end=9910,err=9920) Line
  ReadStatus = ' O.K. after reading ORBA input.'
  call UpCase(Line)
  if (index(Line,'COMP') /= 0) OutFmt2 = 'COMPACT '
  if (index(Line,'FULL') /= 0) OutFmt2 = 'FULL    '
  if (OutFmt2 == 'DEFAULT ') then
    call WarningMessage(1,'Error in ''OUTA'' command?')
    write(LF,*) ' Input line is:'
    write(LF,*) Line
    write(LF,*) ' Did not understand ''OUTA'' input. Ignored.'
  end if
end if
!---  Process MAXO keyword: Max nr of state-specific orbital files produced
if (KeyMAXO) then
  if (DBG) write(6,*) ' MAXORB command:'
  call SetPos(LUInput,'MAXO',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading MAXO input.'
  read(LUInput,*,end=9910,err=9920) MAXORBOUT
  ReadStatus = ' O.K. after reading MAXO input.'
end if
if (DBG) then
  write(6,*) ' Orbital print levels are'
  write(6,*) '   for energy, PRETHR=',PRETHR
  write(6,*) '   for occup , PROTHR=',PROTHR
  write(6,*) ' Orbital listing flag is: ',OutFmt1
  write(6,*) ' Orbital appearance flag: ',OutFmt2
  write(6,*) ' Max nr of state-specific orbital files is ',MAXORBOUT
end if
! ======================================================================
!---  Process OUTORBITALS command (Which kind of orbitals)
if (KeyOUTO) then
  call SetPos(LUInput,'OUTO',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading OUTO input.'
  read(LUInput,*,end=9910,err=9920) Line
  ReadStatus = ' O.K. after reading OUTO input.'
  call UpCase(Line)
  iOrbTyp = 0
  if (index(Line,'AVER') /= 0) iOrbTyp = 1
  if (index(Line,'CANO') /= 0) iOrbTyp = 2
  if (index(Line,'NATU') /= 0) iOrbTyp = 3
  if (index(Line,'SPIN') /= 0) iOrbTyp = 4
  if (iOrbTyp == 0) then
    write(LF,*) ' The line after keyword ''OUTORBITALS'' is'
    write(LF,*) ' not understood. That line begins:'
    write(LF,'(1x,a60)') line(1:60)
    write(LF,*) ' This input is IGNORED.'
  end if
  if (iOrbTyp == 2) IPT2 = 1
  if ((iOrbTyp == 3) .or. (iOrbTyp == 4)) then
    ReadStatus = ' Failure reading nOrbRoot after OUTO keyword.'
    read(LUInput,*,end=9910,err=9920) nOrbRoot
    ReadStatus = ' O.K. after reading nOrbRoot after OUTO keyword.'
  end if
end if
if (DBG) then
  write(6,*) ' OUTORBITALS command specified orbital type ',iOrbTyp
  if (iOrbTyp == 1) write(6,*) ' Meaning: Average'
  if (iOrbTyp == 2) write(6,*) ' Meaning: Canonical'
  if (iOrbTyp == 3) write(6,*) ' Meaning: Natural'
  if (iOrbTyp == 4) write(6,*) ' Meaning: Spin orbitals.'
  if ((iOrbTyp == 3) .or. (iOrbTyp == 4)) write(6,*) ' Max state for printing this orbital type ',nOrbRoot
end if
!---  Process ORDER command (SVC Feb 06)
if (KeyORDE) then
  if (DBG) write(6,*) ' ORDER command was used.'
  call SetPos(LUInput,'ORDE',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure after reading ORDER keyword.'
  read(LUInput,*,end=9910,err=9920) IFORDE
  ReadStatus = ' O.K. after reading ORDER keyword.'
  IORDEM = 1
  if (DBG) write(6,*) ' IFORDE, IORDEM=',IFORDE,IORDEM
end if
!---  Process PRSP command
if (KeyPRSP) then
  if (DBG) write(6,*) ' PRSPIN command was used.'
  ISPDEN = 1
end if
!---  Process IVO command
if (KeyIVO) then
  if (DBG) write(6,*) ' IVO command was used.'
  kIvo = .true.
end if
! ========================================================================
!  If ORBONLY keyword was used, then the JOBIPH file should be used
! only to produce orbital files, then the program stops.
!---  Process ORBO command -( new! generate orbitals from jobiph, only)
if (KeyORBO) then
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
if (KeyALTE) then
  if (DBG) write(6,*) ' ALTER command has been used.'
  call SetPos(LUInput,'ALTE',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure after reading ALTER keyword.'
  read(LUInput,*,end=9910,err=9920) NAlter
  ReadStatus = ' O.K. after reading ALTER keyword.'
  if (NAlter > MaxAlter) then
    write(MaxLab,*) MaxAlter
    MaxLab = adjustl(MaxLab)
    write(LF,*)
    call WarningMessage(2,'Alter: too many orbital pairs.')
    write(LF,*) ' ************* ERROR **************'
    write(LF,*) ' ALTEr: Too many pairs of orbitals '
    write(LF,*) ' to exchange (max '//trim(MaxLab)//').'
    write(LF,*) ' **********************************'
    call Abend()
  end if
  do iAlter=1,NAlter
    ReadStatus = ' Failure reading data after ALTER keyword.'
    read(LUInput,*,end=9910,err=9920) (MAlter(iAlter,i),i=1,3)
    ReadStatus = ' O.K. after reading data after ALTER keyword.'
  end do
  ! (SVC) get absolute orbital values for the alterations so that
  ! iMAlter is symmetry independent
  if (DBG) write(6,*) ' ''Absolute'' iMAlter indices:'
  do iAlter=1,NAlter
    iEnd = 0
    iStart = 1
    do iSym=1,MAlter(iAlter,1)
      iStart = iEnd+1
      iEnd = iEnd+nBas(iSym)
    end do
    iMAlter(iAlter,1) = MAlter(iAlter,2)+iStart-1
    iMAlter(iAlter,2) = MAlter(iAlter,3)+iStart-1
    if (DBG) write(6,'(1x,2I5)') iMAlter(iAlter,1),iMAlter(iAlter,2)
  end do
end if
!---  Process ATOM command (P A Malmqvist Apr 07)
if (KeyATOM) then
  PURIFY = 'ATOM    '
  ISUPSM = 1
  call SetPos(LUInput,'ATOM',Line,iRc)
  call ChkIfKey()
end if
!---  Process LINEAR command (P A Malmqvist Apr 05)
if (KeyLINE) then
  PURIFY = 'LINEAR'
  ISUPSM = 1
  call SetPos(LUInput,'LINE',Line,iRc)
  call ChkIfKey()
end if
if (DBG) write(6,*) ' Purify=',PURIFY

!---  process KSDF command
if (DBG) write(6,*) ' Check if FUNCtional was requested.'
if (KeyFUNC) then
  if (DBG) write(6,*) ' FUNC command was given.'
  DFTFOCK = 'CAS '
  call SetPos(LUInput,'FUNC',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  read(LUInput,*,end=9910,err=9920) Line
  call UpCase(Line)
  if (Line(1:4) == 'ROKS') DFTFOCK = 'ROKS'
  if (Line(1:6) == 'CASDFT') DFTFOCK = 'DIFF'
  read(LUInput,*,end=9910,err=9920) Line
  KSDFT = Line(1:80)
  call UpCase(KSDFT)
  l_casdft = (KSDFT(1:2) == 'T:') .or. (KSDFT(1:3) == 'FT:')
  if (.not. l_casdft) goto 9920
  if ((IPRLOC(1) >= DEBUG) .and. l_casdft) write(6,*) ' MCPDFT with functional:',KSDFT
  ExFac = Get_ExFac(KSDFT)
  !---  Process DFCF command
  if (KeyDFCF) then
    if (DBG) write(6,*) ' DFCF (dens. func. coeff) command was used.'
    call SetPos(LUInput,'DFCF',Line,iRc)
    if (iRc /= _RC_ALL_IS_WELL_) goto 9810
    ReadStatus = ' Failure reading DF coeff after DFCF keyword.'
    read(LUInput,*,end=9910,err=9920) CoefX,CoefR
    ReadStatus = ' O.K. after reading DF coeff after DFCF keyword.'
    if (DBG) then
      write(6,*) ' Density functional exchange coeff, CoefX=',CoefX
      write(6,*) ' Density functional correlation coeff, CoefRE=',CoefR
    end if
    call ChkIfKey()
  end if

  call ChkIfKey()
end if
!---  Process CION command
if (DBG) write(6,*) ' Check if CIONLY case.'
if (KeyCION) then
  if (DBG) write(6,*) ' CIONLY keyword was used.'
  iCIonly = 1
  call SetPos(LUInput,'CION',Line,iRc)
  call ChkIfKey()
end if
!---  Process ROST command
if (DBG) write(6,*) ' Check if ROtSTate case.'
if (KeyROST) then
  if (DBG) write(6,*) ' ROtSTate keyword was used.'
  iRotPsi = 1
  call SetPos(LUInput,'ROST',Line,iRc)
  call ChkIfKey()
end if
!---  Process XMSI command
if (DBG) write(6,*) ' Check if XMSI case.'
if (KeyXMSI) then
  if (DBG) write(6,*) ' XMSI keyword was used.'
  iRotPsi = 1
  IXMSP = 1
  call SetPos(LUInput,'XMSI',Line,iRc)
  call ChkIfKey()
end if
!---  Process CMSI command
if (DBG) write(6,*) ' Check if CMSI case.'
if (KeyCMSI) then
  if (DBG) write(6,*) ' CMSI keyword was used.'
  iRotPsi = 1
  ICMSP = 1
  call SetPos(LUInput,'CMSI',Line,iRc)
  call ChkIfKey()
end if
!---  Process CMSS command
CMSStartMat = 'XMS'
if (KeyCMSS .and. (iCMSP == 1)) then
  if (DBG) write(6,*) ' Reading CMS initial rotation matrix'
  call SetPos(LUInput,'CMSS',Line,iRc)
  Line = Get_Ln(LUInput)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  call ChkIfKey()
  if (DBG) then
    write(6,*) ' Reading CMS starting rotation matrix from'
    write(6,*) trim(Line)
  end if
  if (.not. (Line == 'XMS')) then
    CMSGuessFile = trim(Line)
    CMSStartMat = CMSGuessFile
    call F_Inquire(trim(CMSStartMat),lExists)
    if (.not. lExists) then
      write(LF,'(6X,A,A)') trim(CMSStartMat),' is not found. Use XMS intermediate states as initial guess.'
      CMSStartMat = 'XMS'
    end if
    !call fileorb(Line,CMSStartMat)
  end if
end if
!---  Process CMSO command
if (KeyCMSO .and. (iCMSP == 1)) then
  if (DBG) write(6,*) 'Inputting CMS optimization option'
  call SetPos(LUInput,'CMSO',Line,iRc)
  Line = Get_Ln(LUInput)
  call Upcase(Line)
  if (Line(1:4) == 'NEWT') then
    iCMSOpt = 1
  else if (Line(1:4) == 'JACO') then
    iCMSOpt = 2
  else
    ReadStatus = 'Wrong value assigned to keyword CMSO'
    goto 9920
  end if
  CMSGiveOpt = .true.
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  if (DBG) write(6,*) ' CMS Optimization Option',iCMSOpt
  call ChkIfKey()
end if
!---  Process CMMA command
if (KeyCMMA) then
  if (DBG) write(6,*) ' CMS Max Cylces keyword was given.'
  call SetPos(LUInput,'CMMA',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading data following CMMA keyword.'
  read(LUInput,*,end=9910,err=9920) ICMSIterMax
  ReadStatus = ' O.K. reading data following CMMA keyword.'
  if (DBG) write(6,*) ' Max nr of CMS cycles',ICMSIterMax
  call ChkIfKey()
end if
!---  Process CMMI command
if (KeyCMMI) then
  if (DBG) write(6,*) ' CMS Min Cylces keyword was given.'
  call SetPos(LUInput,'CMMI',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading data following CMMI keyword.'
  read(LUInput,*,end=9910,err=9920) ICMSIterMin
  ReadStatus = ' O.K. reading data following CMMI keyword.'
  if (DBG) write(6,*) ' Min nr of CMS cycles',ICMSIterMin
  call ChkIfKey()
end if
!---  Process CMTH command
if (KeyCMTH) then
  if (DBG) write(6,*) ' CMS Threshold keyword was given.'
  call SetPos(LUInput,'CMTH',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading data following CMTH keyword.'
  read(LUInput,*,end=9910,err=9920) CMSThreshold
  ReadStatus = ' O.K. reading data following CMTH keyword.'
  if (DBG) write(6,*) ' CMS threshold',CMSThreshold
  call ChkIfKey()
end if
!---  Process RFPE command
if (KeyRFPE) then
  if (DBG) write(6,*) ' RFPERT (Response Field Perturbation)'
  RFpert = .true.
  call SetPos(LUInput,'RFPE',Line,iRc)
  call ChkIfKey()
end if
!---  Process NONE command
if (KeyNONE) then
  if (DBG) write(6,*) ' Non-equilibrium response'
  NonEq = .true.
  call SetPos(LUInput,'NONE',Line,iRc)
  call ChkIfKey()
end if
!---  Process RFRO command
if (KeyRFRO) then
  call SetPos(LUInput,'RFRO',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  Line = Get_Ln(LUInput)
  Line(80:80) = '0'
  JPCMROOT = 0
  iall = 0
  ReadStatus = ' Failure reading IPCMROOT after RFROOT keyword.'
  if (nToken(Line) >= 3) then
    read(Line,*,end=9910,err=9920) IPCMROOT,JPCMROOT,iall
    call DWSol_DWRO(LuInput,IPCMROOT,iall)
  else
    read(Line,*,end=9910,err=9920) IPCMROOT
  end if
  ReadStatus = ' O.K. reading IPCMROOT after RFROOT keyword.'

  ! Check that the root value is not changed explicitly by input
  jPCMRoot = iPCMRoot
  call Qpg_iScalar('RF0CASSCF root',Exist)
  if (Exist) then
    call Get_iScalar('RF0CASSCF root',jPCMRoot)
  else
    call Put_iScalar('RF0CASSCF root',iPCMRoot)
  end if

  if (jPCMRoot /= iPCMRoot) then
    !write(6,*) 'iPCMRoot changed by explicitly by input.'

    ! Value changed explicitly by input. Accept the new value.

    call Put_iScalar('RF CASSCF root',iPCMRoot)
    call Put_iScalar('RF0CASSCF root',iPCMRoot)
    call Qpg_dArray('RF CASSCF Vector',Exist,mConf)
    if (Exist) then
      call mma_allocate(RF,mConf,Label='RF')
      RF(:) = 0.0d0
      call Put_dArray('RF CASSCF Vector',RF,mConf)
      call mma_deallocate(RF)
    end if

  else

    ! Value not explicitly changed by input. See if it is
    ! changed by the RunFile, if it exists there.

    call Qpg_iScalar('RF CASSCF root',Exist)
    if (Exist) then
      call Get_iScalar('RF CASSCF root',iPCMRoot)
    else
      call Put_iScalar('RF CASSCF root',iPCMRoot)
    end if

  end if

  if (DBG) then
    write(6,*) ' RFROOT command was given.'
    write(6,*) ' Response field for root number ',IPCMROOT
  end if
  call ChkIfKey()
end if
!---  Process CIRF command
if (KeyCIRF) then
  call SetPos(LUInput,'CIRF',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' O.K. reading after keyword CIRF.'
  read(LUInput,*,end=9910,err=9920) ICIRFROOT
  ReadStatus = ' O.K. reading after keyword CIRF.'
  call ChkIfKey()
  if (DBG) then
    write(6,*) ' CIRFROOT command was given.'
    write(6,*) ' Response field will follow CISE root: ',ICIRFROOT
  end if
end if
!-----------------------------------------------------------------------
if (KeyMCM7) then
# ifndef _HDF5_
  call WarningMessage(2,'MCM7 is given in the input, please make sure to compile Molcas with HDF5 support.')
  goto 9930
# endif
  MCM7 = .true.
  DoNECI = .true.  ! needed to initialise FCIQMC
  totalwalkers = 20000
  RDMsampling%start = 20000
  RDMsampling%n_samples = 10000
  RDMsampling%step = 100
  if (DBG) write(6,*) 'M7 CASSCF activated.'
  if (DBG) write(6,*) 'Decoupled mode not implemented.'
  if (DBG) write(6,*) 'Ignore automatically generated FciInp!'
end if
!-----------------------------------------------------------------------
if (KeyNDPT) then
  tNonDiagStochPT2 = .true.
  IPT2 = 1     ! flag for SXCTL
  if (KeySUPS) then
    write(6,*) 'SUPSymmetry incompatible with NDPT.'
    call Abend()
  end if
  if (KeyPPT2) then
    write(6,*) 'Non-diagonal PT2 incompatible with PPT2.'
    call Abend()
  end if
  if (DBG) write(6,*) 'stoch.-PT2 will be prepared in the current basis.'
  if (DBG) write(6,*) 'Act. Space Fock matrix will be dumped.'
end if
!-----------------------------------------------------------------------
if (KeyPPT2) then
  tPrepStochCASPT2 = .true.
  iOrbTyp = 2  ! pseudo-canonical orbitals
  IPT2 = 1     ! flag for SXCTL
  if (KeySUPS) then
    write(6,*) 'SUPSymmetry incompatible with PPT2.'
    call Abend()
  end if
  if (KeyNDPT) then
    write(6,*) 'Non-diagonal PT2 incompatible with PPT2.'
    call Abend()
  end if
  if (DBG) write(6,*) 'Transforming final orbitals into pseudo-canonical.'
  if (DBG) write(6,*) 'Act. Space Fock matrix will be dumped.'
end if
!-----------------------------------------------------------------------
if (KeyRGRA) then
  tRootGrad = .true.
  if (DBG) write(6,*) 'Orbital gradient for each root is printed.'
end if
!---  Process SSCR command
if (KeySSCR) then
  if (DBG) write(6,*) ' SSCR command was given.'
  call setpos(luinput,'SSCR',line,irc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  line = get_ln(luinput)
  line(80:80) = '0'
  ReadStatus = ' Failure reading after KeySSCR keyword.'
  read(line,*,err=9920,end=9920) norbs,same_orbs
  ReadStatus = ' O.K reading after KeySSCR keyword.'

  if (norbs >= mxOrb) then
    write(6,'(a)',advance='no') 'SSCR error:'
    write(6,*) 'number of spatial orbitals exceeds maximum'
    write(6,'(a,i4)') 'norbs = ',norbs
    write(6,'(a)') new_line('a')
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
      write(6,'(a)',advance='no') 'SSCR error:'
      write(6,*) 'numbers of spatial orbitals do not match'
      write(6,*) 'orb_range_p has length ',size(orb_range_p)
      write(6,*) 'orb_range_q has length ',size(orb_range_q)
      write(6,'(a)') new_line('a')
      call abend()
    end if

    do i=1,norbs
      do j=1,norbs
        if (i < j) then
          if (orb_range_p(i) == orb_range_p(j)) then
            write(6,'(a)',advance='no') 'SSCR error:'
            write(6,*) 'first range contains duplicates.'
            write(6,'(*(i4))') orb_range_p
            write(6,'(a)') new_line('a')
            call abend()
          end if
          if (orb_range_q(i) == orb_range_q(j)) then
            write(6,'(a)',advance='no') 'SSCR error:'
            write(6,*) 'second range contains duplicates.'
            write(6,'(*(i4))') orb_range_q
            write(6,'(a)') new_line('a')
            call abend()
          end if
        end if
      end do
    end do
  else
    do i=1,norbs
      orb_range_p(i) = i
      orb_range_q(i) = i
    end do
  end if
  call ChkIfKey()
end if
!---  Process STAV command
if (KeySTAV .and. KeyCIRO) then
  call WarningMessage(1,'STAVERAGE and CIROOT are incompatible.;The STAVERAGE command will be ignored.')
  KeySTAV = .false.
end if
if (KeySTAV) then
  if (DBG) write(6,*) ' STAVERAGE command was given.'
  call SetPos(LUInput,'STAV',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  Line = Get_Ln(LUInput)
  ReadStatus = ' Failure reading spin after STAVERAGE keyword.'
  read(Line,*,err=9920) NROOTS
  if (NROOTS > MXROOT) then
    write(6,*) 'Error: number of roots exceeds maximum'
    write(6,*) 'NROOTS = ',NROOTS
    write(6,*) 'MXROOT = ',MXROOT
    call AbEnd()
  end if
  ReadStatus = ' O.K. reading spin after STAVERAGE keyword.'
  LROOTS = NROOTS
  do i=1,NROOTS
    iroot(i) = i
    WEIGHT(i) = 1.d0/dble(NROOTS)
  end do
  if (DBG) then
    write(6,*) ' Nr of roots in CI: LROOTS=',LROOTS
    write(6,*) ' Nr of roots optimized by super-CI: NROOTS=',NROOTS
    write(6,*) ' (Equal-weighted)'
  end if
  call ChkIfKey()
end if
!---  Process CIRO command
if (DBG) write(6,*) ' Check for CIROOTS command.'
if (KeyCIRO) then
  if (DBG) write(6,*) ' CIROOTS command was given.'
  call SetPos(LUInput,'CIRO',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  Line = Get_Ln(LUInput)
  !BOR.. Modification 001011
  Line(80:80) = '0'
  ReadStatus = ' Failure reading after CIROOTS keyword.'
  read(Line,*,err=9920,end=9920) NROOTS,LROOTS,iall
  ReadStatus = ' O.K reading after CIROOTS keyword.'
  if (NROOTS > MXROOT) then
    write(6,*) 'Error: number of roots exceeds maximum'
    write(6,*) 'NROOTS = ',NROOTS
    write(6,*) 'MXROOT = ',MXROOT
    call AbEnd()
  end if
  if (iall == 1) then
    do i=1,NROOTS
      iroot(i) = i
      WEIGHT(i) = 1.d0/dble(NROOTS)
    end do
  else
    !BOR.. End modification 001011
    Line = Get_Ln(LUInput)
    ReadStatus = ' Failure reading after CIROOTS keyword.'
    read(Line,*,err=9920,end=9920) (IROOT(I),I=1,NROOTS)
    ReadStatus = ' O.K.after CIROOTS keyword.'
    call dCopy_(mxRoot,[0.0d0],0,WEIGHT,1)
    if (NROOTS == 1) then
      WEIGHT(1) = 1.0d0
    else
      call mma_allocate(Temp1,NROOTS,Label='Temp1')
      Line = Get_Ln(LUInput)
      ReadStatus = ' Failure reading after CIROOTS keyword.'
      read(Line,*,err=9920) (Temp1(i),i=1,NROOTS)
      ReadStatus = ' O.K.after CIROOTS keyword.'
      iSum = 0
      do i=1,NROOTS
        iSum = iSum+Temp1(i)
      end do
      do i=1,NROOTS
        WEIGHT(i) = dble(Temp1(i))/dble(iSum)
      end do
      call mma_deallocate(Temp1)
    end if
  end if
  if (DBG) then
    write(6,*) ' Nr of roots in CI: LROOTS=',LROOTS
    write(6,*) ' Nr of roots optimized by super-CI: NROOTS=',NROOTS
    if (iAll == 1) then
      write(6,*) ' (Equal-weighted)'
    else
      write(6,*) ' Weights:'
      do i1=1,NROOTS,10
        i2 = min(NROOTS,i1+9)
        write(6,'(1x,10f8.4)') (WEIGHT(i),i=i1,i2)
      end do
    end if
  end if
  call ChkIfKey()
end if
!---  Process RLXR command
if (KeyRLXR) then
  call SetPos(LUInput,'RLXR',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' O.K. reading after keyword RLXR.'
  read(LUInput,*,end=9910,err=9920) IRLXROOT
  ReadStatus = ' O.K. reading after keyword RLXR.'
  call ChkIfKey()
  if (DBG) then
    write(6,*) ' RLXROOT command was given.'
    write(6,*) ' State for SLAPAF to handle: ',IRLXROOT
  end if
  if (.not. any(iRoot(1:LROOTS) == IRLXROOT)) then
    write(6,*) ' The selected root is not among those available.'
    call AbEnd()
  end if
end if

!IgorS 29-4-2010 Begin
!---  Process MDRL command
if (KeyMDRL) then
  if (KeyRLXR) then
    write(6,*) ' RLXROOT keyword was given before MDRLXROOT.'
    write(6,*) ' Since these keywords are mutually exclusive'
    write(6,*) ' please check the input and read the manual.'
    goto 9910
  end if
  call SetPos(LUInput,'MDRL',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  call Qpg_iScalar('Relax CASSCF root',RlxRCheck)
  if (RlxRCheck) then
    call Get_iScalar('Relax CASSCF root',iRlxRoot)
    if (DBG) then
      write(6,*) ' An existing relax root was found.'
      write(6,*) ' The MDRLXROOT value is ignored.'
    end if
  else
    ReadStatus = ' Failure reading relaxroot number IRLXROOT.'
    read(LUInput,*,end=9910,err=9920) IRLXROOT
    ReadStatus = ' O.K. after reading relaxroot number IRLXROOT.'
    call ChkIfKey()
    if (.not. any(iRoot(1:LROOTS) == IRLXROOT)) then
      write(6,*) ' The selected root is not among those available.'
      call AbEnd()
    end if
  end if
  if (DBG) then
    write(6,*) ' MDRLxroot command was given.'
    write(6,*) ' DYNAMIX will follow the root: ',IRLXROOT
  end if
end if
!IgorS End

!---  Process CISE command
if (DBG) write(6,*) ' Check for CISELECT command.'
if (KeyCISE) then
  if (DBG) then
    write(6,*) ' CISELECT keyword was given.'
    write(6,*) ' This input is awkward. Let''s find up'
    write(6,*) ' a better way to do things.'
  end if
  ICICH = 1
  call SetPos(LUInput,'CISE',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  do i=1,NROOTS
    ReadStatus = ' Failure reading after CISELECT keyword.'
    read(LUInput,*,end=9910,err=9920) kRef
    ReadStatus = ' O.K. reading after CISELECT keyword.'
    if (kRef > mxRef) then
      call WarningMessage(1,'CISElect input is wrong.')
      write(LF,*) 'Number of CSF''s in CiSelect is out of bounds'
      write(LF,'(a,i3,a,i3)') 'Specified:',kRef,', Max is',mxRef
      write(LF,'(a,i3)') 'Standard fixup, value set to',mxRef
      kRef = mxRef
    end if
    ReadStatus = ' Failure reading after CISELECT keyword.'
    read(LUInput,*,end=9910,err=9920) (ICI(i,iRef),iRef=1,kRef)
    read(LUInput,*,end=9910,err=9920) (CCI(i,iRef),iRef=1,kRef)
    ReadStatus = ' O.K. reading after CISELECT keyword.'
    dSum = 0.0d0
    do iRef=1,kRef
      dSum = dSum+CCI(i,iRef)**2
    end do
    do iRef=1,kRef
      CCI(i,iRef) = CCI(i,iRef)/sqrt(dSum)
    end do
    do iRef=kRef+1,mxRef
      CCI(i,iRef) = 0.0d0
      ICI(i,iRef) = 0
    end do
  end do
  call ChkIfKey()
end if

!---  Process ALPH command
if (KeyALPH) then
  if (DBG) write(6,*) ' The ALPH keyword was used.'
  call SetPos(LUInput,'ALPH',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading data after ALPH keyword.'
  read(LUInput,*,end=9910,err=9920) iAlphaBeta
  ReadStatus = ' O.K. after reading data after ALPH keyword.'
  if (iAlphaBeta < 0) iAlphaBeta = -1
  if (iAlphaBeta > 0) iAlphaBeta = 1
  if (DBG) then
    if (iAlphaBeta == 1) write(6,*) ' Read alpha orbitals from UHF'
    if (iAlphaBeta == -1) write(6,*) ' Read beta orbitals from UHF'
  end if
  call ChkIfKey()
end if

! =========   Input source for orbitals: ===============================
! INVEC=0 is used to indicate if any source of orbitals has been
! identified.
! iOrbData=0 is used to indicate if any source of orbital type
! information -- inactive, ras1, etc -- has been identified.

if (DBG) then
  write(6,*) ' Check for input source of orbitals.'
  write(6,*) 'KeyCORE,KeyJOBI,KeyLUMO:',KeyCORE,KeyJOBI,KeyLUMO
end if
! Handle multiple keywords:
if (KeyCORE) then
  KeyJOBI = .false.
  KeyLUMO = .false.
else if (KeyJOBI) then
  KeyLUMO = .false.
end if
! Only one of these should have been selected.
if (DBG) write(6,*) 'KeyCORE,KeyJOBI,KeyLUMO:',KeyCORE,KeyJOBI,KeyLUMO

! CORE is probably becoming obsolete.
!---  Process CORE command
if (KeyCORE) then
  if (DBG) write(6,*) ' CORE command was used.'
  if (IPRLEV >= VERBOSE) write(LF,*) ' Start orbitals will be computed from scratch.'
  INVEC = 1
end if

!---  Process JOBI command
if (KeyJOBI) then
  if (DBG) write(6,*) ' JOBIPH command was used.'
  call f_Inquire('JOBOLD',lExists)
  if (lExists) then
    if (IPRLEV >= VERBOSE) write(LF,*) ' Orbitals will be taken from the old jobiph file.'
    INVEC = 3
  else
    call f_Inquire(IPHNAME,lExists)
    if (lExists) then
      if (IPRLEV >= VERBOSE) write(LF,*) ' Orbitals will be taken from the jobiph file.'
      INVEC = 3
    end if
  end if
  if (INVEC == 0) then
    if (IPRLEV >= TERSE) then
      call WarningMessage(2,'JOBIPH input is wrong.')
      write(6,*) ' Keyword JOBIPH was used, but the ''JOBOLD'' file'
      write(6,*) ' does not exist. The ''JOBIPH'' file named'
      write(6,*) IPHNAME
      write(6,*) 'also does not exist. This is a fatal error.'
      goto 9930
    end if
  end if
end if

!---  Process H5OR command
if (KeyH5OR) then
# ifdef _HDF5_
  KeyLUMO = .false.
  KeyTYPE = .false.
  iOverwr = merge(1,0,any([KeyRAS1,KeyRAS2,KeyRAS3,KeyFROZ,KeyINAC,KeyDELE]))
  mh5id = mh5_open_file_r(StartOrbFile)
  ! read basic attributes
  call mh5_fetch_attr(mh5id,'NSYM',NSYM_L)
  if (nsym /= nsym_l) then
    write(6,*) 'Number of symmetries on HDF5 file does not'
    write(6,*) 'match the number of symmetries on the'
    write(6,*) 'RunFile, calculation will stop now.'
    call Quit(_RC_INPUT_ERROR_)
  end if
  call mh5_fetch_attr(mh5id,'NBAS',NBAS_L)
  ierr = 0
  do isym=1,nsym
    if (nbas(isym) /= nbas_l(isym)) ierr = 1
  end do
  if (ierr == 1) then
    write(6,*) 'Number of basis functions on HDF5 file does not'
    write(6,*) 'match the number of basis functions on the'
    write(6,*) 'RunFile, calculation will stop now.'
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
  if (KeyH5CI) then
    if (mh5_exists_dset(mh5id,'CI_VECTORS') .and. mh5_exists_attr(mh5id,'NROOTS')) then
      call mh5_fetch_attr(mh5id,'NROOTS',lroots_l)
      if (lroots_l /= lroots) then
        write(6,*) 'Number of CI roots on file does not'
        write(6,*) 'match the number requested by the user,'
        write(6,*) 'so no CI vectors will be read from'
        write(6,*) StartOrbFile
        iCIRST = 0
      else
        iCIRST = 1
      end if
    else
      write(6,*) 'The required fields CI_VECTORS and/or NROOTS'
      write(6,*) 'are missing from the HDF5 file supplied by'
      write(6,*) 'the user. As a result, to continue,'
      write(6,*) 'no CI vectors will be read from'
      write(6,*) StartOrbFile
      iCIRST = 0
    end if
  end if
  call mh5_close_file(mh5id)
# else
  write(6,*) 'The format of the start orbital file was'
  write(6,*) 'specified by the user as HDF5, but this'
  write(6,*) 'is not implemented in this installation.'
  call Quit(_RC_INPUT_ERROR_)
# endif
end if

!---  Process LUMO command
InfoLbl = '        '
if (KeyLUMO) then
  if (DBG) then
    write(6,*) ' LUMORB command was used.'
    write(6,*) ' Name of orbital file, StartOrbFile='
    write(6,*) StartOrbFile
  end if
  call ChkVec(StartOrbFile,inporb_version,NSYM_L,NBAS_L,NORB_L,InfoLbl,iRc1)
  if (iRc1 /= _RC_ALL_IS_WELL_) then
    if (IPRLEV >= TERSE) then
      call WarningMessage(1,'LUMORB input error.')
      write(6,*) ' Keyword LUMORB used with file name StartOrbFile='
      write(6,*) StartOrbFile
      write(6,*) ' but that file cannot be used. Perhaps it does'
      write(6,*) ' not exist?'
    end if
    iOverWr = 1
  else
    if (DBG) then
      write(6,*) ' The file may be used to read input orbitals.'
      write(6,*) ' It is of type INPORB version ',inporb_version
      write(6,*) ' The information label is: ',InfoLbl
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
      write(6,*) ' ERROR: Start orbital file name is '
      write(6,*) StartOrbFile
      write(6,*) ' That file is a valid orbital file.'
      write(6,*) ' Version:',inporb_version
      write(6,*) ' But some information does not match.'
      if (IERR == 1) then
        write(6,*) ' In the file, nr of symmetries is =',NSYM_L
        write(6,*) ' but according to the runfile, it is=',NSYM
      else if (IERR == 2) then
        write(6,*) ' In the file, nr of basis functions/symm is'
        write(6,'(1x,8I5)') (NBAS_L(I),I=1,NSYM)
        write(6,*) ' but RUNFILE claims it is'
        write(6,'(1x,8I5)') (NBAS(I),I=1,NSYM)
      end if
      write(6,*) ' Is it an old file left in workspace by mistake?'
      goto 9930
    end if

    ! This also implies that information on orbital types could be
    ! taken from typeindex on orbital file:
    if (('I' .in. to_upper(trim(InfoLbl))) .and. (.not. any([KeyRAS1,KeyRAS2,KeyRAS3,KeyFROZ,KeyINAC,KeyDELE]))) then
      iOrbData = 3
      iOverWr = 0
      if (DBG) then
        write(6,*) ' This means we may take orbital specifications'
        write(6,*) ' from the file, so set iOrbData=3, iOverWr=0.'
        write(6,*) ' The orbital spaces are read from typeindex.'
      end if
      ! We will also take the opportunity to find the orbital spaces size
      ! according to typeindex, for possible need below:
      call mma_allocate(type,mxOrb,Label='Type')
      LuStartOrb = 19
      call RdVec(StartOrbFile,LuStartOrb,'IA',NSYM_L,NBAS_L,NBAS_L,Dummy,Dummy,Dummy,type,myTitle,0,iErr)
      call tpidx2orb(NSYM_L,NBAS_L,type,NFRO_L,NISH_L,NRS1_L,NRS2_L,NRS3_L,NSSH_L,NDEL_L)
      call mma_deallocate(type)
      if (DBG) then
        write(6,*) ' From RDTPIDX, we get:'
        write(6,'(1x,A16,8I4)') ' NBAS_L:',(NBAS_L(I),I=1,NSYM_L)
        write(6,'(1x,A16,8I4)') ' NORB_L:',(NORB_L(I),I=1,NSYM_L)
        write(6,'(1x,A16,8I4)') ' NFRO_L:',(NFRO_L(I),I=1,NSYM_L)
        write(6,'(1x,A16,8I4)') ' NISH_L:',(NISH_L(I),I=1,NSYM_L)
        write(6,'(1x,A16,8I4)') ' NRS1_L:',(NRS1_L(I),I=1,NSYM_L)
        write(6,'(1x,A16,8I4)') ' NRS2_L:',(NRS2_L(I),I=1,NSYM_L)
        write(6,'(1x,A16,8I4)') ' NRS3_L:',(NRS3_L(I),I=1,NSYM_L)
        write(6,'(1x,A16,8I4)') ' NSSH_L:',(NSSH_L(I),I=1,NSYM_L)
        write(6,'(1x,A16,8I4)') ' NDEL_L:',(NDEL_L(I),I=1,NSYM_L)
      end if
    else
      iOverWr = 1
    end if
  end if
  INVEC = 2
end if
if (DBG) write(6,*) ' The INVEC    code is now',INVEC
if (DBG) write(6,*) ' The iOrbData code is now',iOrbData

! ======================================================================

!---  Process TYPEindex command
if (DBG) write(6,*) ' Was TYPEINDEX requested?'
if (KeyTYPE) then
  call SetPos(LUInput,'TYPE',Line,iRc)
  call ChkIfKey()
  if (DBG) then
    write(6,*) ' TYPEINDEX command was used.'
    write(6,*) ' The size of orbital spaces should be read from'
    write(6,*) ' typeindex in starting orbital file.'
  end if
  if ((index(InfoLbl,'i') > 0) .or. (index(InfoLbl,'I') > 0)) then
    iOrbData = 3
    iOverwr = 0
    if (IPRLEV >= VERBOSE) write(LF,*) ' Orbital specification will be taken from orbital file'
    call mma_allocate(type,mxOrb,Label='Type')
    LuStartOrb = 19
    call RdVec(StartOrbFile,LuStartOrb,'IA',NSYM_L,NBAS_L,NBAS_L,Dummy,Dummy,Dummy,type,myTitle,0,iErr)
    call tpidx2orb(NSYM_L,NBAS_L,type,NFRO_L,NISH_L,NRS1_L,NRS2_L,NRS3_L,NSSH_L,NDEL_L)
    call mma_deallocate(type)
    IERR = 0
    if (NSYM_L /= NSYM) IERR = 1
    if (IERR == 0) then
      do ISYM=1,NSYM
        if (NBAS(ISYM) /= NBAS_L(ISYM)) IERR = 1
      end do
    end if
    if ((IERR /= 0) .or. DBG) then
      write(LF,*)
      write(LF,'(6X,A)') 'Specifications read from runfile:'
      write(LF,'(6X,A)') '----------------------------------------'
      write(LF,*)
      write(LF,'(6X,A,T47,8I4)') 'Symmetry species',(iSym,iSym=1,NSYM)
      write(LF,'(6X,A,T47,8I4)') 'Number of basis functions',(NBAS(iSym),iSym=1,NSYM)
      write(LF,*)
      write(LF,'(6X,A)') 'Specifications read from orbital file:'
      write(LF,'(6X,A)') '----------------------------------------'
      write(LF,*)
      write(LF,'(6X,A,T47,8I4)') 'Symmetry species',(iSym,iSym=1,NSYM_L)
      write(LF,'(6X,A,T47,8I4)') 'Frozen orbitals',(NFRO_L(iSym),iSym=1,NSYM_L)
      write(LF,'(6X,A,T47,8I4)') 'Inactive orbitals',(NISH_L(iSym),iSym=1,NSYM_L)
      write(LF,'(6X,A,T47,8I4)') 'RAS1 orbitals',(NRS1_L(iSym),iSym=1,NSYM_L)
      write(LF,'(6X,A,T47,8I4)') 'RAS2 orbitals',(NRS2_L(iSym),iSym=1,NSYM_L)
      write(LF,'(6X,A,T47,8I4)') 'RAS3 orbitals',(NRS3_L(iSym),iSym=1,NSYM_L)
      write(LF,'(6X,A,T47,8I4)') 'Secondary orbitals',(NSSH_L(iSym),iSym=1,NSYM_L)
      write(LF,'(6X,A,T47,8I4)') 'Deleted orbitals',(NDEL_L(iSym),iSym=1,NSYM_L)
      write(LF,'(6X,A,T47,8I4)') 'Number of basis functions',(NBAS_L(iSym),iSym=1,NSYM_L)
      write(LF,*)
    end if
    if ((IERR /= 0) .and. (IPRLEV >= TERSE)) then
      write(6,*) ' Orbital specifications were to be read from'
      write(6,*) ' orbital file, but there is mismatch with'
      write(6,*) ' some data on the runfile!'
      write(6,*) ' Orbital file name is:'
      write(6,*) StartOrbFile
      iOrbData = 0
    else
      do ISYM=1,NSYM
        NFRO(ISYM) = NFRO_L(ISYM)
        NISH(ISYM) = NISH_L(ISYM)
        NRS1(ISYM) = NRS1_L(ISYM)
        NRS2(ISYM) = NRS2_L(ISYM)
        NRS3(ISYM) = NRS3_L(ISYM)
        NSSH(ISYM) = NSSH_L(ISYM)
        NDEL(ISYM) = NDEL_L(ISYM)
      end do
    end if
  else
    if (DBG) write(6,*) ' Typeindex cannot be read!'
  end if
end if
if (DBG) write(6,*) ' The iOrbData code is now',iOrbData

! ======================================================================
! Explicit orbital sizes input:
! Save a copy of current iorbdata first:
iod_save = iorbdata
!---  Process FROZ command
if (KeyFROZ) then
  if (DBG) write(6,*) ' FROZEN keyword was given.'
  call SetPos(LUInput,'FROZ',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading after FROZEN keyword.'
  read(LUInput,*,end=9910,err=9920) (NFRO(ISYM),ISYM=1,NSYM)
  ReadStatus = ' O.K. reading after FROZEN keyword.'
  call Get_iScalar('nSym',i)
  call Put_iArray('nFro',nFro,i)
  if (DBG) then
    write(6,*) ' Nr of Frozen orbitals requested:'
    write(6,'(1x,8i5)') (NFRO(i),i=1,NSYM)
  end if
  IORBDATA = 1
end if
!---  Process INAC command
if (KeyINAC) then
  if (DBG) write(6,*) ' INACTIVE keyword was given.'
  call SetPos(LUInput,'INAC',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading after INACTIVE keyword.'
  read(LUInput,*,end=9910,err=9920) (NISH(ISYM),ISYM=1,NSYM)
  ReadStatus = ' O.K. reading after INACTIVE keyword.'
  if (DBG) then
    write(6,*) ' Nr of Inactive orbitals requested:'
    write(6,'(1x,8i5)') (NISH(i),i=1,NSYM)
  end if
  IORBDATA = 1
end if

!---  Process RAS1 command
if (KeyRAS1) then
  if (DBG) write(6,*) ' RAS1 keyword was given.'
  call SetPos(LUInput,'RAS1',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading after RAS1 keyword.'
  read(LUInput,*,end=9910,err=9920) (NRS1(ISYM),ISYM=1,NSYM)
  ReadStatus = ' O.K. reading after RAS1 keyword.'
  if (DBG) then
    write(6,*) ' Nr of RAS1 orbitals requested:'
    write(6,'(1x,8i5)') (NRS1(i),i=1,NSYM)
  end if
  IORBDATA = 1
end if

!---  Process RAS2 command
if (KeyRAS2) then
  if (DBG) write(6,*) ' RAS2 keyword was given.'
  call SetPos(LUInput,'RAS2',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading after RAS2 keyword.'
  read(LUInput,*,end=9910,err=9920) (NRS2(ISYM),ISYM=1,NSYM)
  ReadStatus = ' O.K. reading after RAS2 keyword.'
  if (DBG) then
    write(6,*) ' Nr of Ras2 orbitals requested:'
    write(6,'(1x,8i5)') (NRS2(i),i=1,NSYM)
  end if
  IORBDATA = 1
end if

!---  Process RAS3 command
if (KeyRAS3) then
  if (DBG) write(6,*) ' RAS3 keyword was given.'
  call SetPos(LUInput,'RAS3',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading after RAS3 keyword.'
  read(LUInput,*,end=9910,err=9920) (NRS3(ISYM),ISYM=1,NSYM)
  ReadStatus = ' O.K. reading after RAS3 keyword.'
  if (DBG) then
    write(6,*) ' Nr of Ras3 orbitals requested:'
    write(6,'(1x,8i5)') (NRS3(i),i=1,NSYM)
  end if
  IORBDATA = 1
end if

!---  Process GASS command
if (KeyGASS) then
  if (DBG) write(6,*) 'GAS is actived'
  call setpos(luinput,'GASS',line,irc)
  if (irc /= _RC_ALL_IS_WELL_) goto 9810
  read(luinput,*,end=9910,err=9920) NGAS
  do igas=1,ngas
    read(luinput,*,end=9910,err=9920) (ngssh(igas,isym),isym=1,nsym)
    read(luinput,*,end=9910,err=9920) (igsoccx(igas,mm),mm=1,2)
  end do
  iDoGas = .true.
  iorbdata = 1
end if

!---  Process DELE command
if (KeyDELE) then
  if (DBG) write(6,*) ' DELETED keyword was given.'
  call SetPos(LUInput,'DELE',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading after DELETED keyword.'
  read(LUInput,*,end=9910,err=9920) (NDEL(ISYM),ISYM=1,NSYM)
  ReadStatus = ' O.K. after reading after DELETED keyword.'
  if (DBG) then
    write(6,*) ' Nr of Deleted orbitals requested:'
    write(6,'(1x,8i5)') (NDEL(i),i=1,NSYM)
  end if
  IORBDATA = 1
end if

if ((IORBDATA == 1) .and. (IPRLEV >= VERBOSE)) write(LF,*) ' Orbital specification was read from input.'
if ((IOD_SAVE == 3) .and. (IORBDATA == 1)) then
  ! See if the input matches the values on file
  IERR = 0
  do ISYM=1,NSYM
    if (NFRO(ISYM) /= NFRO_L(ISYM)) IERR = 1
    if (NISH(ISYM) /= NISH_L(ISYM)) IERR = 1
    if (NRS1(ISYM) /= NRS1_L(ISYM)) IERR = 1
    if (NRS2(ISYM) /= NRS2_L(ISYM)) IERR = 1
    if (NRS3(ISYM) /= NRS3_L(ISYM)) IERR = 1
    if (NSSH(ISYM) /= NSSH_L(ISYM)) IERR = 1
    if (NDEL(ISYM) /= NDEL_L(ISYM)) IERR = 1
  end do
  if (IERR == 0) then
    if ((IORBDATA == 1) .and. (IPRLEV >= VERBOSE)) then
      write(LF,*) ' However, input matches the typeindex on the'
      write(LF,*) ' starting orbitals file. Therefore, accept'
      write(LF,*) ' the typeindex information for sorting.'
      iOrbData = 3
      iOverWr = 0
    end if
  end if
end if
if (IORBDATA == 3) then
  do ISYM=1,NSYM
    NFRO(ISYM) = NFRO_L(ISYM)
    NISH(ISYM) = NISH_L(ISYM)
    NRS1(ISYM) = NRS1_L(ISYM)
    NRS2(ISYM) = NRS2_L(ISYM)
    NRS3(ISYM) = NRS3_L(ISYM)
    NSSH(ISYM) = NSSH_L(ISYM)
    NDEL(ISYM) = NDEL_L(ISYM)
  end do
end if
! ======================================================================
! If IORBDATA is still 0, lets hope there is information on the runfile.
! Exception: If this is a CIRESTART, it must be taken from the JOBIPH
! (or JOBOLD) file.
if (IORBDATA == 0) then
  if (IPRLEV >= VERBOSE) write(LF,*) ' No explicit orbital specs in user input.'
  if (KeyCIRE) then
    if (IPRLEV >= VERBOSE) then
      write(LF,*) ' This is a CIRESTART case, so take them from'
      write(LF,*) ' the JOBIPH or JOBOLD file.'
    end if
    IORBDATA = 2
    if (IPRLEV >= VERBOSE) write(LF,*) ' Orbital specs taken from JOBIPH or JOBOLD.'
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
    call WR_RASSCF_Info(JobOld,2,iAd19,NACTEL,ISPIN,NSYM,STSYM,NFRO,NISH,NASH,NDEL,NBAS,mxSym,lJobH1,(LenIn+8)*mxOrb,NCONF,lJobH2, &
                        2*72,JobTit,4*18*mxTit,POTNUCDUMMY,LROOTS,NROOTS,IROOT,mxRoot,NRS1,NRS2,NRS3,NHOLE1,NELEC3,IPT2,WEIGHT)
  else
    if (IPRLEV >= VERBOSE) then
      write(LF,*) ' This is not a CIRESTART case, so take them from'
      write(LF,*) ' the RUNFILE.'
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
NTOT = 0
NTOT1 = 0
NTOT2 = 0
NO2M = 0
NISHT = 0
NASHT = 0
NDELT = 0
NFROT = 0
NSEC = 0
NORBT = 0
NTOT3 = 0
NTOTSP = 0
NTOT4 = 0
NRS1T = 0 ! for RASSCF
NRS2T = 0
NRS3T = 0
!NGSSH_tot(:) = Zero
!do igas=1,ngas
!  NGSSH_tot(igas) = SUM(NGSSH(IGAS,1:NSYM))
!end do
!if (dbg) then
!  write(6,*) 'NGSSH_tot(igas):'
!  write(6,*) (NGSSH_tot(igas),igas=1,ngas)
!end if
do ISYM=1,NSYM
  NTOT = NTOT+NBAS(ISYM)
  NTOT1 = NTOT1+NBAS(ISYM)*(NBAS(ISYM)+1)/2
  NTOT2 = NTOT2+NBAS(ISYM)**2
  NO2M = max(NO2M,NBAS(ISYM)**2)
  NRS1T = NRS1T+NRS1(ISYM)  ! for RAS
  NRS2T = NRS2T+NRS2(ISYM)
  NRS3T = NRS3T+NRS3(ISYM)
  NFROT = NFROT+NFRO(ISYM)
  NISHT = NISHT+NISH(ISYM)
  NASHT = NASHT+NASH(ISYM)
  NDELT = NDELT+NDEL(ISYM)
  NSEC = NSEC+NSSH(ISYM)
  NORBT = NORBT+NORB(ISYM)
  NTOT3 = NTOT3+(NORB(ISYM)+NORB(ISYM)**2)/2
  NTOTSP = NTOTSP+(NASH(ISYM)*(NASH(ISYM)+1)/2)
  NTOT4 = NTOT4+NORB(ISYM)**2
end do
NACPAR = (NASHT+NASHT**2)/2
NACPR2 = (NACPAR+NACPAR**2)/2
! NASHT is called NAC in some places:
NAC = NASHT
! Same, NISHT, NIN:
NIN = NISHT
NFR = NFROT
if (DBG) write(6,*) ' The iOrbData code is now',iOrbData
! ======================================================================
! Compute effective nuclear charge.
! Identical to nr of protons for conventional basis sets only, not ECP.
call Get_iScalar('Unique atoms',nNuc)
call mma_allocate(ENC,nNuc,Label='ENC')
call Get_dArray('Effective nuclear Charge',ENC,nNuc)
TEffNChrg = 0.0d0
call mma_allocate(Stab,nNuc,Label='Stab')
call Get_iArray('nStab',Stab,nNuc)
do i=1,nNuc
  TEffNChrg = TEffNChrg+ENC(i)*dble(nSym/Stab(i))
end do
call mma_deallocate(Stab)
call mma_deallocate(ENC)
if (DBG) write(6,*) ' Effective nuclear charge is TEffNChrg=',TEffNChrg
TotChrg = 0.0d0
if (DBG) write(6,*) ' Set TotChrg=',TotChrg
!---  Process NACT command
if (KeyNACT) then
  ! Cannot set the number of inactive orbitals if explicitly given or
  ! or if there is symmetry
  if (KeyCHAR .and. (KeyINAC .or. (NSYM > 1))) then
    if (IPRLEV >= TERSE) call WarningMessage(1,'CHARGE and NACTEL are only compatible if INACTIVE is not given and the '// &
                                             'symmetry group is C1. Hence the CHARGE command will be ignored.')
    KeyCHAR = .false.
  end if
  if (DBG) write(6,*) ' NACTEL keyword was given.'
  call SetPos(LUInput,'NACT',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  Line = Get_Ln(LUInput)
  ReadStatus = ' Failure reading data after NACTEL keyword.'
  read(Line,*,err=9920,end=123) NACTEL,NHOLE1,NELEC3
  goto 124
123 read(Line,*,err=9920,end=9920) NACTEL
  NHOLE1 = 0
  NELEC3 = 0
124 continue
  ReadStatus = ' O.K. after reading data after NACTEL keyword.'
  if (DBG) write(6,*) ' NACTEL,NHOLE1,NELEC3:',NACTEL,NHOLE1,NELEC3
  ! Only set total charge here if not explicitly given
  if (.not. KeyCHAR) then
    TotChrg = TEffNChrg-dble(2.0d0*(NISHT+NFROT)+NACTEL)
    if (DBG) write(6,*) ' TotChrg=',TotChrg
  end if
  call Put_iScalar('nActel',NACTEL)
  call ChkIfKey()
end if
!---  Process CHARGE command
if (KeyCHAR) then
  if (DBG) write(6,*) ' CHARGE command was given.'
  call SetPos(LUInput,'CHAR',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  Line = Get_Ln(LUInput)
  ReadStatus = ' Failure reading charge after CHARGE keyword.'
  read(Line,*,err=9920) NCHRG
  ReadStatus = ' O.K. reading charge after CHARGE keyword.'
  if (DBG) write(6,*) ' Total charge is ',NCHRG
  TotChrg = dble(NCHRG)
  call ChkIfKey()
  ! If both CHAR and NACT where given, set the inactive and secondary
  ! orbitals accordingly
  if (KeyNACT) then
    NISHT_old = NISHT
    NISHT = (int(TEffNChrg-TotChrg+0.5d0)-NACTEL)/2-NFROT
    NISH(1) = NISHT
    NIN = NISHT
    NSEC = NSEC-(NISHT-NISHT_old)
    NSSH(1) = NSEC
  end if
end if
! The NINT function is unreliable on Cygwin gfortran, use INT:
! Nr of electrons should be positive integer, so this is probably safe:
NACTEL = int(TEffNChrg-TotChrg+0.5d0)-2*(NISHT+NFROT)
call Put_iScalar('nActel',NACTEL)
if (DBG) then
  write(6,*) ' Compute NActEl from  other data:'
  write(6,*) '     TEffNChrg=',TEffNChrg
  write(6,*) '       TotChrg=',TotChrg
  write(6,*) '       NISHT  =',NISHT
  write(6,*) '       NFROT  =',NFROT
  write(6,*) ' Resulting NActEl=',NActEl
end if
!---  Process RASSCF command
if (KeyRASS) then
  if (DBG) write(6,*) ' RASSCF keyword was given.'
  call SetPos(LUInput,'RASS',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  Line = Get_Ln(LUInput)
  ReadStatus = ' Failure reading data after RASSCF keyword.'
  read(Line,*,err=9920,end=9920) NHOLE1,NELEC3
  ReadStatus = ' O.K. reading data after RASSCF keyword.'
  call ChkIfKey()
end if
! ======================================================================
!---  Process SPIN command
if (DBG) write(6,*) ' Determine spin value:'
if (KeySPIN) then
  if (DBG) write(6,*) ' SPIN command was given.'
  call SetPos(LUInput,'SPIN',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  Line = Get_Ln(LUInput)
  ReadStatus = ' Failure reading spin after SPIN keyword.'
  read(Line,*,err=9920) ISPIN
  ReadStatus = ' O.K. reading spin after SPIN keyword.'
  if (DBG) then
    write(6,*) ' Spin multiplicity is ',ISPIN
    S = 0.5d0*dble(ISPIN-1)
    write(6,*) '  i.e. SPIN=',S
  end if
  call ChkIfKey()
else
  ! If ISPIN has not been set, use some value that may have been set
  ! in runfile (e.g., UHF, or previous calculation)
  if (DBG) write(6,*) ' No SPIN command was given.'
  call qpg_iscalar('ISPIN',lExists)
  if (lExists) then
    call get_iscalar('ISPIN',ISPIN)
    if (DBG) write(6,*) ' Runfile has ISPIN=',ISPIN
  else
    if (DBG) write(6,*) ' Runfile does not know ISPIN.'
    call qpg_dscalar('UHFSPIN',lExists)
    if (lExists) then
      call get_dscalar('UHFSPIN',SUHF)
      ISPIN = nint(1.0d0+2.0d0*SUHF)
      if (DBG) write(6,*) ' Runfile has UHFSPIN=',SUHF
    else
      ! or, last chance fallback, guess on singlet or doublet:
      if (DBG) write(6,*) ' Runfile does not know UHFSPIN.'
      ISPIN = 1
      if (NACTEL /= 2*(NACTEL/2)) ISPIN = 2
      if (IPRLEV >= TERSE) then
        call WarningMessage(1,'Had to guess the spin.')
        write(6,*) ' Warning: no input and no reliable source'
        write(6,*) ' for the spin multiplicity.'
        write(6,*) ' Guess ISPIN=',ISPIN
      end if
    end if
  end if
end if
call put_iscalar('ISPIN',ISPIN)
! If spin is zero, do not compute and print spin density:
if (ISPIN == 1) ISPDEN = 0
! ======================================================================
if (KeyDMPO) DumpOnly = .true.
! ======================================================================
if (KeyReOr) then
  if (DBG) write(6,*) 'Orbital Reordering (REOR) is activated'
  call setpos(luinput,'REOR',line,irc)
  if (irc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading ReOrFlag after REOR keyword.'
  read(luinput,*,end=9910,err=9920) ReOrFlag
  ReadStatus = ' O.K. reading ReOrFlag after REOR keyword.'
  if ((ReOrFlag < -1) .or. (ReOrFlag == 1)) then
    call WarningMessage(2,"Invalid flag for reordering. n==0: Don't reorder. n>=2: User defined permutation with n changed "// &
                        'elements. n==-1: Use GAS sorting scheme. ')
    goto 9930
  else if (ReOrFlag >= 2) then
    call mma_allocate(ReOrInp,ReOrFlag)
    ReadStatus = ' Failure reading ReOrInp after REOR keyword.'
    read(luinput,*,end=9910,err=9920) (ReOrInp(i),i=1,ReOrFlag)
    ReadStatus = ' O.K. reading ReOrInp after REOR keyword.'
  else if ((ReOrFlag == -1) .and. (.not. KeyGASS)) then
    call WarningMessage(2,'If GAS is not used, a permutation for orbital reordering has to be specified.')
    goto 9930
  end if
end if
if (KeyORTH) then
  call setpos(luinput,'ORTH',line,irc)
  if (irc /= _RC_ALL_IS_WELL_) goto 9810
  read(luinput,*,end=9910,err=9920) ON_scheme_inp
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
    goto 9930
  end if
end if
if (KeyPERI) write_orb_per_iter = .true.
!---  Process NECI commands
if (KeyNECI) then
  if (DBG) write(6,*) 'NECI is actived'
  DoNECI = .true.

  if (KeyDMPO) then
    call WarningMessage(2,'NECI and DMPOnly are mutually exclusive.')
    goto 9930
  end if
  !---------------------------------------------------------------------
  if (KeyEMBD) then
    DoEmbdNECI = .true.
#ifndef _NECI_
    call WarningMessage(2,'EmbdNECI is given in input, so the embedded NECI should be used. Unfortunately MOLCAS was not '// &
                        'compiled with embedded NECI. Please use -DNECI=ON for compiling or use an external NECI.')
#endif
  end if
  !---------------------------------------------------------------------
  if (KeyWRMA) then
    WRMA = .true.
    if (DBG) write(6,*) 'DMAT/PSMAT/PAMAT will be dumped.'
  end if
  !---------------------------------------------------------------------
  if (KeyGUGA) then
    tGUGA_in = .true.
    if (DBG) write(6,*) 'spin-free GUGA-NECI RDMs are actived'
    if (.not. KeyNECI) then
      call WarningMessage(2,'GUGA requires NECI keyword!')
      goto 9930
    end if
  end if
  !--- This block is to process the DEFINEDET
  if (KeyDEFI) then
    call setpos(luinput,'DEFI',line,irc)
    call mma_allocate(buffer,2000)
    if (irc /= _RC_ALL_IS_WELL_) goto 9810
    ReadStatus = ' Failure reading Definedet.'
    read(luinput,'(A)',end=9910,err=9920) buffer
    ReadStatus = ' O.K. reading Definedet.'
    call mma_allocate(definedet,len_trim(buffer))
    definedet(:) = trim(buffer)
    call mma_deallocate(buffer)
    write(6,*) 'definedet read in proc_inp of size:',nactel
    write(6,*) definedet
  end if
  if (KeyTOTA) then
    call setpos(luinput,'TOTA',line,irc)
    if (irc /= _RC_ALL_IS_WELL_) goto 9810
    ReadStatus = ' Failure reading totalwalkers.'
    read(luinput,*,end=9910,err=9920) totalwalkers
    ReadStatus = ' O.K. reading totalwalkers.'
  else
    call WarningMessage(2,'TOTAlwalkers required for NECI.')
    goto 9930
  end if
  if (count([KeyRDML,(KeyRDMS .and. KeyCALC)]) /= 1) then
    call WarningMessage(2,'RDMLinspace, and (RDMSamplingiters + CALCrdmonfly) are mutually exclusive, but one is required.')
    goto 9930
  else if (KeyRDML) then
    call setpos(luinput,'RDML',line,irc)
    if (irc /= _RC_ALL_IS_WELL_) goto 9810
    read(luinput,*,end=9910,err=9920) RDMsampling%start,RDMsampling%n_samples,RDMsampling%step
  else if (KeyRDMS .or. KeyCALC) then
    if (.not. (KeyRDMS .and. KeyCALC)) then
      call WarningMessage(2,'RDMSamplingiters and CALCrdmonfly are both required.')
      goto 9930
    end if
    call setpos(luinput,'RDMS',line,irc)
    if (irc /= _RC_ALL_IS_WELL_) goto 9810
    read(luinput,*,end=9910,err=9920) start

    call setpos(luinput,'CALC',line,irc)
    if (irc /= _RC_ALL_IS_WELL_) goto 9810
    read(luinput,*,end=9910,err=9920) length,step
    RDMsampling = t_RDMsampling(start,length/step,step)
  end if
  if (KeyDIAG) then
    call setpos(luinput,'DIAG',line,irc)
    if (irc /= _RC_ALL_IS_WELL_) goto 9810
    read(luinput,*,end=9910,err=9920) diagshift
  end if
  if (KeyTIME) then
    call setpos(luinput,'TIME',line,irc)
    if (irc /= _RC_ALL_IS_WELL_) goto 9810
    read(luinput,*,end=9910,err=9920) Time
  end if
  if (KeyNMCY) then
    call setpos(luinput,'NMCY',line,irc)
    if (irc /= _RC_ALL_IS_WELL_) goto 9810
    read(luinput,*,end=9910,err=9920) nmCyc
  end if
  if (KeyREAL) then
    call setpos(luinput,'REAL',line,irc)
    if (irc /= _RC_ALL_IS_WELL_) goto 9810
    read(luinput,*,end=9910,err=9920) realspawncutoff
  end if
  if (KeyTRIA) then
    call setpos(luinput,'TRIA',line,irc)
    if (irc /= _RC_ALL_IS_WELL_) goto 9810
    read(luinput,*,end=9910,err=9920) trial_wavefunction
  end if
  if (KeyPOPS) then
    call setpos(luinput,'POPS',line,irc)
    if (irc /= _RC_ALL_IS_WELL_) goto 9810
    read(luinput,*,end=9910,err=9920) pops_trial
  end if
  if (KeySEMI) then
    call setpos(luinput,'SEMI',line,irc)
    if (irc /= _RC_ALL_IS_WELL_) goto 9810
    read(luinput,*,end=9910,err=9920) semi_stochastic
  end if
  if (KeyMEMO) then
    call setpos(luinput,'MEMO',line,irc)
    if (irc /= _RC_ALL_IS_WELL_) goto 9810
    read(luinput,*,end=9910,err=9920) memoryfacspawn
  end if
  !call fciqmc_option_check(iDoGas, nGSSH, iGSOCCX)
end if

if (KeyCCCI) then
  if (DBG) write(6,*) 'CC-CI is actived'
  Do_CC_CI = .true.

  if (KeyDMPO) then
    call WarningMessage(2,'CC-CI and DMPOnly are mutually exclusive.')
    goto 9930
  end if
end if

! ======================================================================
!---  Process SYMM command
if (KeySYMM) then
  if (DBG) write(6,*) ' SYMMETRY command was given.'
  call SetPos(LUInput,'SYMM',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  Line = Get_Ln(LUInput)
  ReadStatus = ' Failure reading symmetry index after SYMM keyword.'
  read(Line,*,err=9920) STSYM
  ReadStatus = ' O.K. reading symmetry index after SYMM keyword.'
  if (DBG) write(6,*) ' State symmetry index ',STSYM
  call ChkIfKey()
else
  ! If STSYM has not been set, normally it should be defaulted to 1.
  ! Exception: if this is a high-spin OS case, these often require STSYM /= 1:
  STSYM = 1
  if (ISPIN == NASHT+1) then
    do ISYM=1,NSYM
      NA = NASH(ISYM)
      if (NA /= 2*(NA/2)) STSYM = MUL(STSYM,ISYM)
    end do
  end if
end if
call put_iscalar('STSYM',STSYM)
if (DBG) write(6,*) ' State symmetry STSYM=',STSYM

! ======================================================================

!---  Process CIRE command
if (KeyCIRE) then
  if (DBG) write(6,*) ' CIRESTART keyword was given.'
  ICIRST = 1
end if

!---  Process HOME command (root homing in SXCI part)
if (KeyHOME) then
  SXSEL = 'HOMING  '
  if (DBG) write(6,*) ' HOME (Root Homing) keyword was given.'
  call SetPos(LUInput,'HOME',Line,iRc)
  call ChkIfKey()
end if

!---  Process SUPS command
if (KeySUPS) then
  if (DBG) write(6,*) ' SUPS (Supersymmetry) keyword was given.'
  call SetPos(LUInput,'SUPS',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  call mma_allocate(Temp1,mxOrb,Label='Temp1')
  ISUPSM = 1
  iOffset = 0
  do iSym=1,nSym
    ReadStatus = ' Failure reading data following SUPS keyword.'
    read(LUInput,*,end=9910,err=9920) nGrp
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
if (KEYHEXS) then
  if (DBG) write(6,*) ' HEXS (Highly excited states) keyword was given. '
  call SetPos(LUInput,'HEXS',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  if ((I_ELIMINATE_GAS_MOLCAS /= 0) .and. (I_ELIMINATE_GAS_MOLCAS /= 2)) then
    call WarningMessage(2,'HEXS keyword defined more than once')
    goto 9810
  end if
  I_ELIMINATE_GAS_MOLCAS = I_ELIMINATE_GAS_MOLCAS+1
  ReadStatus = ' Failure reading data following HEXS keyword.'
  read(LUInput,*,end=9910,err=9920) N_ELIMINATED_GAS_MOLCAS
  ReadStatus = ' O.K. after reading data following HEXS keyword.'
  ReadStatus = ' Failure reading data following HEXS keyword.'
  read(LUInput,*,end=9910,err=9920) (IELIMINATED_IN_GAS_MOLCAS(I),I=1,N_ELIMINATED_GAS_MOLCAS)
  ReadStatus = ' O.K. after reading data following HEXS keyword.'
end if

! --- Process DEXS command
if (KEYDEXS) then
  if (DBG) write(6,*) ' DEXS (Doubly excited states) keyword was given. '
  call SetPos(LUInput,'DEXS',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  if (I_ELIMINATE_GAS_MOLCAS > 1) then
    call WarningMessage(2,'DEXS keyword defined more than once')
    goto 9810
  end if
  I_ELIMINATE_GAS_MOLCAS = I_ELIMINATE_GAS_MOLCAS+2
  ReadStatus = ' Failure reading data following DEXS keyword.'
  read(LUInput,*,end=9910,err=9920) N_2ELIMINATED_GAS_MOLCAS
  ReadStatus = ' O.K. after reading data following DEXS keyword.'
  ReadStatus = ' Failure reading data following DEXS keyword.'
  read(LUInput,*,end=9910,err=9920) (I2ELIMINATED_IN_GAS_MOLCAS(I),I=1,N_2ELIMINATED_GAS_MOLCAS)
  ReadStatus = ' O.K. after reading data following DEXS keyword.'
end if

!---  Process HROO command ---
if (KEYHROO) then
  if (DBG) write(6,*) ' HROO (Hidden roots) keyword was given. '
  call SetPos(LUInput,'HROO',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading data following HROO keyword.'
  read(LUInput,*,end=9910,err=9920) hRoots
  ReadStatus = ' O.K. after reading data following HROO keyword.'
end if
!---  Process NKEE command ---
if (KEYNKEE) then
  if (DBG) write(6,*) ' NKEE (nr of kept vectors) keyword was given. '
  call SetPos(LUInput,'NKEE',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading data following NKEE keyword.'
  read(LUInput,*,end=9910,err=9920) n_keep
  ReadStatus = ' O.K. after reading data following NKEE keyword.'
  if (n_keep < lRoots) then
    call WarningMessage(2,'nkeep must be at least equal to the number of roots')
    call Quit(_RC_INPUT_ERROR_)
  end if
end if

!---  Process CLEA command ---
if (KeyCLEA) then
  if (DBG) write(6,*) ' CLEAN (Orbital Cleaning) keyword.'
  if (DBG) write(6,*) ' (Awkward input -- replace??).'
  call SetPos(LUInput,'CLEA',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  call mma_allocate(Temp1,mxOrb,Label='Temp1')
  call mma_allocate(Temp2,mxOrb,Label='Temp2')
  call mma_allocate(Temp3,mxOrb,Label='Temp3')
  nClean = 0
  do iSym=1,nSym
    nClean = nClean+nBas(iSym)**2
  end do
  call mma_allocate(Cleanmask,nClean,Label='CleanMask')
  iOffset = 0
  do iSym=1,nSym
    mBas = nBas(iSym)
    do i=1,mBas
      ii = (i-1)*mBas
      do j=1,mBas
        ij = j+ii+iOffset
        CleanMask(ij) = 0
      end do
    end do
    ReadStatus = ' Failure reading data following CLEAN keyword.'
    read(LUInput,*,end=9910,err=9920) nGrp
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
if (KeyCHOL) then
  if (DBG) write(6,*) ' CHOLESKY keyword was given.'
  DoCholesky = .true.
  call SetPos(LUInput,'CHOL',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  call Cho_rasscf_rdInp(.true.,LUInput)
end if

!---  Process CHOI command (Cholesky Custom Input, F.Aquilante)
! Cholesky with user-defined settings.
if (KeyCHOI) then
  if (DBG) write(6,*) ' CHOINPUT keyword was given.'
  DoCholesky = .true.
  call SetPos(LUInput,'CHOI',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  call Cho_rasscf_rdInp(.false.,LUInput)
end if

!---  Process ITER command
if (KeyITER) then
  if (DBG) write(6,*) ' ITERATIONS keyword was given.'
  call SetPos(LUInput,'ITER',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading data following ITER keyword.'
  read(LUInput,*,end=9910,err=9920) MAXIT,ITMAX
  ReadStatus = ' O.K. reading data following ITER keyword.'
  if (DBG) then
    write(6,*) ' Max nr of RASSCF (macro) iterations MAXIT=',MAXIT
    write(6,*) ' Max nr of orbital optimization iter ITMAX=',ITMAX
    if (KeyCION) write(6,*) ' (Irrelevant, since the CIONLY keyword was also given)'
  end if
  call ChkIfKey()
end if

!---  Process CRPR command
if (KeyCRPR) then
  if (DBG) write(6,*) ' CRPR (Core Projector) keyword was used.'
  call SetPos(LUInput,'CRPR',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading core orbital after CRPR keyword.'
  read(LUInput,*,end=9910,err=9920) ITCORE,CORESHIFT
  ReadStatus = ' O.K. reading core orbital after CRPR keyword.'
  IfCRPR = .true.
  if (DBG) write(6,*) ' Core orbital is number ITCORE'
  call ChkIfKey()
end if

!---  Process LEVS command
if (KeyLEVS) then
  if (DBG) write(6,*) ' LEVS (Level Shift) keyword was used.'
  call SetPos(LUInput,'LEVS',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading level shift after LEVSHFT keyword.'
  read(LUInput,*,end=9910,err=9920) LVSHFT
  ReadStatus = ' O.K. reading level shift after LEVSHFT keyword.'
  if (DBG) write(6,*) ' Level shift LVSHFT (Re*8!!) =',LVSHFT
  call ChkIfKey()
end if

!---  Process THRS command
if (KeyTHRS) then
  if (DBG) write(6,*) ' THRS (Thresholds) command was used.'
  call SetPos(LUInput,'THRS',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading thresholds after THRS keyword.'
  read(LUInput,*,end=9910,err=9920) THRE,THRTE,THRSX
  ReadStatus = ' O.K. after reading thresholds after THRS keyword.'
  !tbp, may 2013: no altering of thre with Cholesky
  !tbp if (DoCholesky .and. (IPRLEV >= TERSE)) then
  !tbp   write(6,*) '*** Detected Cholesky or RI/DF calculation'
  !tbp   write(6,*) ' BUT user specified value of ThrE will be used. ThrE= ',THRE
  !tbp end if
  if (DBG) then
    write(6,*) ' Threshold for energy change,    THRE =',THRE
    write(6,*) ' Threshold for orbital rotation, THRTE=',THRTE
    write(6,*) ' Threshold for max BLB element,  THRSX=',THRSX
  end if
  call ChkIfKey()
end if

!---  Process TIGH command
if (KeyTIGH) then
  if (DBG) write(6,*) ' TIGHT (Tight CI convergence)  used.'
  KTIGHT = 1
  call SetPos(LUInput,'TIGH',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading data after TIGHT keyword.'
  read(LUInput,*,end=9910,err=9920) THREN,THFACT
  ReadStatus = ' O.K. reading data after TIGHT keyword.'
  if (DBG) then
    write(6,*) ' CI energy threshold in 1st iter, THREN=',THREN
    write(6,*) ' CI Threshold/Energy change =    THFACT=',THFACT
  end if
  call ChkIfKey()
end if

! Use of quasi-newton ci/orbital coupling?
! Commands QUNE or NOQUNE:
if (KeyNOQU) then
  NQUNE = 0
  if (DBG) write(6,*) ' NOQUNE keyword: QUNE is disabled.'
  call SetPos(LUInput,'NOQU',Line,iRc)
  call ChkIfKey()
else if (KeyQUNE) then
  if (DBG) write(6,*) ' QUNE keyword: QUNE is enabled.'
  NQUNE = 1
  call SetPos(LUInput,'QUNE',Line,iRc)
  call ChkIfKey()
else
  ! Default is to use QUNE, unless this is some kind of DFT:
  if (KeyFUNC) then
    NQUNE = 0
    if (DBG) write(6,*) ' DFT calculation: QUNE is disabled.'
  else
    NQUNE = 1
    if (DBG) write(6,*) ' QUNE is enabled by default.'
  end if
end if

!---  Process CIMX command
if (KeyCIMX) then
  if (DBG) write(6,*) ' Keyword CIMX was used.'
  call SetPos(LUInput,'CIMX',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading data after CIMX keyword.'
  read(LUInput,*,end=9910,err=9920) MAXJT
  ReadStatus = ' O.K. after reading data after CIMX keyword.'
  if (DBG) write(6,*) ' Max nr of CI iterations MAXJT=',MAXJT
  call ChkIfKey()
end if

!---  Process SDAV command
if (KeySDAV) then
  if (DBG) write(6,*) ' SDAV (Size of explicit Hamiltonian matrix)'
  call SetPos(LUInput,'SDAV',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading data after SDAV keyword.'
  read(LUInput,*,end=9910,err=9920) NSEL
  ReadStatus = ' O.K. after reading data after SDAV keyword.'
  if (DBG) write(6,*) ' Size is NSEL=',NSEL
  call ChkIfKey()
end if

!---  Process OFEM commands for Orbital-Free embedding
if (KeyOFEM) then
  if (DBG) write(6,*) ' OFEM (Orbital-Free Embedding activated)'
  Do_OFemb = .true.
  call SetPos(LuInput,'OFEM',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading data after OFEM keyword.'
  read(LUInput,'(A)',end=9910,err=9920) OFE_KSDFT
  ReadStatus = ' O.K. after reading data after OFEM keyword.'
  call UpCase(OFE_KSDFT)
  OFE_KSDFT = adjustl(OFE_KSDFT)
  write(6,*)
  write(6,*) '  --------------------------------------'
  write(6,*) '   Orbital-Free Embedding Calculation'
  write(6,*) '  --------------------------------------'
  if (OFE_KSDFT(1:4) == 'LDTF') then
    write(6,*) '    T_nad potential   : Thomas-Fermi    '
  else
    write(6,*) '    T_nad potential   : ',OFE_KSDFT(1:4)
  end if
  if (KEonly) then
    write(6,*) '    Exc_nad potential :  None           '
  else
    write(6,*) '    Exc_nad potential : ',OFE_KSDFT(6:10)
  end if
  write(6,*) '  --------------------------------------'
  write(6,*)
  DFTFOCK = 'SCF '
end if
if (KeyFTHA) then
  call SetPos(LuInput,'FTHA',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading data after FTHA keyword.'
  read(LUInput,*,end=9910,err=9920) ThrFThaw
  ReadStatus = ' O.K. after reading data after FTHA keyword.'
end if
if (KeyDFMD) then
  call SetPos(LuInput,'DFMD',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading data after DFMD keyword.'
  read(LUInput,*,end=9910,err=9920) dFMD,Xsigma
  ReadStatus = ' O.K. after reading data after DFMD keyword.'
  !write(6,'(A,F6.3)') '  Fract. correl. potent. DFMD=',dFMD
  !write(6,*) '  --------------------------------------'
  if (dFMD+Xsigma < 0.0d0) then
    write(6,*) ' *** Warning: arguments to DFMD must be nonnegative!'
    write(6,*) ' ***          I will take their ABS !!! '
    dFMD = abs(dFMD)
    Xsigma = abs(Xsigma)
  end if
end if

!---  Process BKAP command for BK type of approximation
!     (Giovanni Li Manni J.:GLMJ) Nov 2011
if (KeyBKAP) then
  DoBKAP = .true.
  call SetPos(LUInput,'BKAP',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading data after BKAP keyword.'
  read(LUInput,*,end=9910,err=9920) NGASBK
  read(LUInput,*,end=9910,err=9920) (IOCCPSPC(IGAS,1),IGAS=1,NGASBK)
  read(LUInput,*,end=9910,err=9920) (IOCCPSPC(IGAS,2),IGAS=1,NGASBK)
  ReadStatus = ' O.K. reading data after BKAP keyword.'
  if (DBG) then
    write(6,*) ' BKAP: BK-type of approximation in action'
    write(6,*) ' Min and Max for subspace with exact Hamiltonian '
    write(6,*) ' =============================================== '
    write(6,*) ' NGASBK :',NGASBK
    write(6,*) '              Min. Occ.      Max. Occ.           '
    do IGAS=1,NGASBK
      write(6,'(A,I2,10X,I3,9X,I3)') '   GAS',IGAS,IOCCPSPC(IGAS,1),IOCCPSPC(IGAS,2)
    end do
  end if
end if

!---  Process SPLI command for SplitCAS calculations
!     (Giovanni Li Manni J.:GLMJ)
if (KeySPLI) then
  if (DBG) write(6,*) ' SPLI (Activation SplitCAS)'
  DoSplitCAS = .true.
  EnerSplit = .true.
  iDimBlockA = 0
  !*** The energy gap (GapSpli) is in eV ****
  GapSpli = 13.61d0
  lrootSplit = 1
  thrSplit = 1.0d-6
  MxIterSplit = 50
  call SetPos(LUInput,'SPLI',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
end if

!---  Process NUSP command for Numerical SplitCAS param. (GLMJ)
if (KeyNUSP) then
  if (DBG) write(6,*) ' NUSP - Manual Setting of Numerical SplitCAS Param.'
  EnerSplit = .false.
  NumSplit = .true.
  call SetPos(LUInput,'NUSP',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading data after NUSP keyword.'
  read(LUInput,*,end=9910,err=9920) lrootSplit,iDimBlockA,MxIterSplit
  ReadStatus = ' O.K. reading data after NUSP keyword.'
  !ead(LUInput,*,end=9910,err=9920) lrootSplit
  !ead(LUInput,*,end=9910,err=9920) iDimBlockA
  !ead(LUInput,*,end=9910,err=9920) MxIterSplit
  ReadStatus = ' Failure reading data after NUSP keyword.'
  read(LUInput,*,end=9910,err=9920) ThrSplit
  ReadStatus = ' O.K. reading data after NUSP keyword.'
  if (DBG) then
    write(6,*) ' Root to be opt. in SplitCAS = ',lrootSplit
    write(6,*) ' AA block size in SplitCAS = ',iDimBlockA
    write(6,*) ' Max iteration in SplitCAS = ',MxIterSplit
    write(6,*) ' Root to be opt. in SplitCAS = ',ThrSplit
  end if
end if

!---  Process ENSP command for Energetical SplitCAS param. (GLMJ)
if (KeyENSP) then
  if (DBG) write(6,*) ' ENSP - Manual Setting of Energetical SplitCAS Param.'
  EnerSplit = .true.
  NumSplit = .false.
  iDimBlockA = 0
  call SetPos(LUInput,'ENSP',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading data after ENSP keyword.'
  read(LUInput,*,end=9910,err=9920) lrootSplit,GapSpli,MxIterSplit
  ReadStatus = ' O.K. reading data after ENSP keyword.'
  !read(LUInput,*,end=9910,err=9920) lrootSplit
  !read(LUInput,*,end=9910,err=9920) iDimBlockA
  !read(LUInput,*,end=9910,err=9920) MxIterSplit
  ReadStatus = ' Failure reading data after ENSP keyword.'
  read(LUInput,*,end=9910,err=9920) ThrSplit
  ReadStatus = ' O.K. reading data after ENSP keyword.'
  if (DBG) write(6,*) ' Root to be opt. in SplitCAS = ',lrootSplit
  if (DBG) write(6,*) ' Energy gap in SplitCAS = ',GapSpli
  if (DBG) write(6,*) ' Max iteration in SplitCAS = ',MxIterSplit
  if (DBG) write(6,*) ' Root to be opt. in SplitCAS = ',ThrSplit
end if

!---  Process PESP command for Percentage SplitCAS param. (GLMJ)
if (KeyPESP) then
  if (DBG) write(6,*) ' PESP - Manual Setting of Percentage SplitCAS Param.'
  EnerSplit = .false.
  PerSplit = .true.
  NumSplit = .false.
  iDimBlockA = 0
  call SetPos(LUInput,'PESP',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading data after PESP keyword.'
  read(LUInput,*,end=9910,err=9920) lrootSplit,PercSpli,MxIterSplit
  read(LUInput,*,end=9910,err=9920) ThrSplit
  ReadStatus = ' O.K. reading data after PESP keyword.'
  if (DBG) write(6,*) ' Root to be opt. in SplitCAS = ',lrootSplit
  if (DBG) write(6,*) ' Percentage in SplitCAS = ',PercSpli
  if (DBG) write(6,*) ' Max iteration in SplitCAS = ',MxIterSplit
  if (DBG) write(6,*) ' Root to be opt. in SplitCAS = ',ThrSplit
end if

!---  Process FOSP command for First Order SplitCAS Approx. (GLMJ)
if (KeyFOSP) then
  if (DBG) write(6,*) ' FOSP - First Order SplitCAS Approx.'
  FOrdSplit = .true.
end if

!---  Process OPTO keyword: Optimal Output for RASSCF/CASPT2
!                           optimizations - GG Nov 2008.
if (KeyOPTO) then
  if (DBG) then
    write(6,*) ' OPTO keyword was used.'
    write(6,*) '(Optimal Output for RASSCF/CASPT2)'
  end if
  lOPTO = .true.
  call SetPos(LUInput,'OPTO',Line,iRc)
  call ChkIfKey()
end if

!---  Process SXDAmp command
if (KeySXDA) then
  if (DBG) write(6,*) ' SXDAMPING was requested.'
  call SetPos(LUInput,'SXDA',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading data after SXDAMP keyword.'
  read(LUInput,*,end=9910,err=9920) SXDamp
  ReadStatus = ' O.K. after reading data after SXDAMP keyword.'
  if (DBG) write(6,*) ' Parameter SXDamp=',SXDamp
  call ChkIfKey()
end if

!---  Process LOWM command
if (KeyLOWM) then
  if (DBG) then
    write(6,*) ' LOWM keyword was used to force the CI routines'
    write(6,*) ' to use Slater Determinants for low M and not M=S'
  end if
  LOWMS = 1
  call SetPos(LUInput,'LOWM',Line,iRc)
  call ChkIfKey()
end if

!---  Process LOWD keyword: Turn on Lowdin orthonormalization of CMOs
if (KeyLOWD) then
  if (DBG) then
    write(6,*) ' LOWDIN orthonormalization was requested'
    write(6,*) ' but from Jan 12 2010 that is default anyway.'
  end if
  Lowdin_ON = .true.
  call SetPos(LUInput,'LOWD',Line,iRc)
  call ChkIfKey()
end if

!---  Process PRWF command
if (KeyPRWF) then
  if (DBG) write(6,*) ' The PRWF keyword was used.'
  call SetPos(LUInput,'PRWF',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading data after PRWF keyword.'
  read(LUInput,*,end=9910,err=9920) PRWTHR
  ReadStatus = ' O.K. after reading data after PRWF keyword.'
  if (DBG) write(6,*) ' Print CI coefficients larger than PRWTHR=',PRWTHR
  call ChkIfKey()
end if

!---  Process PRSD command
if (KeyPRSD) then
  if (DBG) write(6,*) ' The PRSD keyword was used.'
  call SetPos(LUInput,'PRSD',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  if (DBG) write(6,*) ' Print determinant expansions of CSFs'
  call ChkIfKey()
end if

!---  Process FCIDUMP command
if (KeyFCID) then
  ! activate the DMRG interface in RASSCF (dummy here since we stop after FCIDUMP)
  DOFCIDUMP = .true.
  if (.not. KeyDMRG) KeyDMRG = .true.
  call SetPos(LUInput,'FCID',Line,iRc)
  call ChkIfKey()
end if

#if ! defined (_ENABLE_BLOCK_DMRG_) && ! defined (_ENABLE_CHEMPS2_DMRG_)
! ======================================================================
!          start of QCMaquis DMRG input section
! ======================================================================
#ifdef _DMRG_
if (keyDMRG .or. doDMRG) then
  if (.not. doDMRG) then
    call SetPos(LUInput,'DMRG',Line,iRc)
    call ChkIfKey()
  end if
  !> DMRG flag
  doDMRG = .true.
  LRras2_dmrg(1:8) = 0
  !> LRras2 = Ras2 as the default
  do i=1,nsym
    LRras2_dmrg(i) = NRS2(i)
  end do
  !> initial guess setup
  guess_dmrg(1:7) = 'DEFAULT'
  call mma_allocate(initial_occ,nrs2t,nroots); initial_occ = 0
end if

!---  Process RGIN command (QCMaquis Custom Input)
if (KeyRGIN) then
  if (DBG) write(6,*) ' RGINPUT keyword was given.'
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  call SetPos(LUInput,'RGIN',Line,iRc)

  if (.not. KeyDMRG) then
    if (.not. doDMRG) then
      call WarningMessage(2,'Error in input processing.')
      write(6,*) ' PROC_INP: the keyword DMRG is not present but'
      write(6,*) ' is required to enable the DMRG internal keyword'
      write(6,*) ' section RGInput.'
      iRc = _RC_INPUT_ERROR_
      Go to 9900
    end if
  end if

  nr_lines = 0
  call qcmaquis_rdinp(LuInput,1,nr_lines)
  call SetPos(LUInput,'RGIN',Line,iRc)
  call qcmaquis_rdinp(LuInput,2,nr_lines)

end if

!---  Process SOCC command (state occupation for initial guess in DMRG)
if (KeySOCC) then
  if (DBG) write(6,*) ' SOCC keyword was given.'
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  call SetPos(LUInput,'SOCC',Line,iRc)

  if (keyDMRG .or. doDMRG) then
    call socc_dmrg_rdinp(LuInput,initial_occ,nrs2t,nroots)
    guess_dmrg(1:7) = 'HF     '
  end if
end if
#endif

!-- Leon: Process NEVP(t2prep) keyword, prepare for 4-RDM calculation
!    for (CD)-DMRG-NEVPT2
if (KeyNEVP) then
  if (DBG) write(6,*) ' NEVP(t2prep) keyword was given.'
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
# ifdef _DMRG_
  if ((.not. KeyDMRG) .and. (.not. doDMRG)) then
    call WarningMessage(2,'Error in input processing.')
    write(6,*) ' PROC_INP: the keyword DMRG is not present or'
    write(6,*) ' DMRG not activated but is required to enable'
    write(6,*) ' the NEVP keyword.'
    iRc = _RC_INPUT_ERROR_
    Go to 9900
  end if

  DoNEVPT2Prep = .true.
  call SetPos(LUInput,'NEVP',Line,iRc)
  Line = Get_Ln(LUInput)
  call UpCase(Line)
  if (index(Line,'EVRD') /= 0) then
    call WarningMessage(2,'Warning,EvRDM keyword is deprecated.')
    write(6,*) 'RDM evaluation is done in NEVPT2 module now'
  end if

# else
  call WarningMessage(2,'Error in input processing.')
  write(6,*) 'MOLCAS was compiled without QCMaquis support.'
  write(6,*) 'Thus, no DMRG-NEVPT2 calculations are possible.'
  iRc = _RC_INPUT_ERROR_
  Go to 9900
# endif
end if

#ifdef _DMRG_
!> sanity checks
!> a. DMRG requested but mandatory keywords not set at all
if ((KeyDMRG .or. doDMRG) .and. ((.not. KeyRGIN) .and. (.not. as_solver_inp_proc))) then
  call WarningMessage(2,'Error in input processing.')
  write(6,*) ' PROC_INP: the keyword RGINput is not present but'
  write(6,*) ' is required for QCMaquis DMRG calculations in order'
  write(6,*) ' to set compulsory DMRG internal parameters as'
  write(6,*) ' for example:'
  write(6,*) ' max_bond_dimension, conv_thresh and nsweeps.'
  write(6,*) ' See the QCMaquis manual for further details.'
  iRc = _RC_INPUT_ERROR_
  Go to 9900
end if
!> b. DMRG requested so check that ALL mandatory keywords have been set
if ((KeyDMRG .or. doDMRG) .and. (KeyRGIN .or. as_solver_inp_proc)) then
  nr_lines = dmrg_input%nr_qcmaquis_input_lines
  call qcmaquis_rdinp(LuInput,3,nr_lines)
  if (nr_lines <= 0) then
    iRc = _RC_INPUT_ERROR_
    Go to 9900
  end if
end if
#endif
! ======================================================================
!          end of QCMaquis DMRG input section
! ======================================================================
#endif

!---  Process FARO command
if (KeyFARO) DoFaro = .true.

!---  Process NOCA command
if (DBG) write(6,*) ' Check if NOCALC case.'
if (KeyNOCA) then
  if (DBG) write(6,*) ' NOCALC keyword was used.'
  INOCALC = 1
  call SetPos(LUInput,'NOCA',Line,iRc)
  call ChkIfKey()
end if

!---  Process SAVE command
if (DBG) write(6,*) ' Check if SAVE_EXP case.'
if (KeySAVE) then
  if (DBG) write(6,*) ' SAVE_EXP keyword was used.'
  ISAVE_EXP = 1
  call SetPos(LUInput,'SAVE',Line,iRc)
  call ChkIfKey()
end if

!---  Process EXPA command
if (DBG) write(6,*) ' Check if EXPAND case.'
if (KeyEXPA) then
  if (DBG) write(6,*) ' EXPAND keyword was used.'
  IEXPAND = 1
  call SetPos(LUInput,'EXPA',Line,iRc)
  call ChkIfKey()
end if

!---  Process DMRG command
#if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_)
if (KeyDMRG) then
  ! NN.14 FIXME: When DMRG option is disabled at compilation,
  !       this should give an error, but just ignored for the time.
  if (DBG) write(6,*) ' DMRG (Use DMRG algorithm instead of FCI)'
  call SetPos(LUInput,'DMRG',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading data after DMRG keyword.'
  read(LUInput,*,end=9910,err=9920) MxDMRG
  ReadStatus = ' O.K. after reading data after DMRG keyword.'
  if (DBG) write(6,*) ' Nr. of states=',MxDMRG
  DoBlockDMRG = .true.
  call ChkIfKey()
end if

!---  Process 3RDM command
if (Key3RDM) then
  if (DBG) write(6,*) ' 3RDM (Compute 3RDM for DMRG-Cu4-CASPT2)'
  Do3RDM = .true.
# ifdef _ENABLE_CHEMPS2_DMRG_
  iOrbTyp = 2
  IPT2 = 1
  write(6,*) 'CHEMPS2> 3-RDM and F4-RDM require PseudoCanonical orbitals'
  write(6,*) 'CHEMPS2> Automatically set: OUTOrbitals = CANOnical'
  if (KeySUPS) then
    write(6,*) 'CHEMPS2> Bug using SYPSym and 3RDM, disable SUPSym'
    call Abend()
  end if
# endif
  call SetPos(LUInput,'3RDM',Line,iRc)
  call ChkIfKey()
end if

#ifdef _ENABLE_CHEMPS2_DMRG_
!---  Process DAVT command
if (KeyDAVT) then
  call SetPos(LUInput,'DAVT',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading data after DAVT keyword.'
  read(LUInput,*,end=9910,err=9920) davidson_tol
  ReadStatus = ' O.K. after reading data after DAVT keyword.'
  call ChkIfKey()
end if

!---  Process CHRE command
if (KeyCHRE) then
  if (DBG) write(6,*) ' Restart in CheMPS2'
  chemps2_restart = .true.
  call SetPos(LUInput,'CHRE',Line,iRc)
  call ChkIfKey()
end if

!---  Process CHBL command
if (KeyCHBL) then
  call SetPos(LUInput,'CHBL',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading data after CHBL keyword.'
  read(LUInput,*,end=9910,err=9920) chemps2_blb
  ReadStatus = ' O.K. after reading data after CHBL keyword.'
  call ChkIfKey()
end if

!---  Process MXSW command
if (KeyMXSW) then
  call SetPos(LUInput,'MXSW',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading data after MXSW keyword.'
  read(LUInput,*,end=9910,err=9920) max_sweep
  ReadStatus = ' O.K. after reading data after MXSW keyword.'
  call ChkIfKey()
end if

!---  Process NOIS command
if (KeyNOIS) then
  call SetPos(LUInput,'NOIS',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading data after NOIS keyword.'
  read(LUInput,*,end=9910,err=9920) chemps2_noise
  ReadStatus = ' O.K. after reading data after NOIS keyword.'
  call ChkIfKey()
end if

!---  Process DMRE command
if (KeyDMRE) then
  call SetPos(LUInput,'DMRE',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading data after DMRE keyword.'
  read(LUInput,*,end=9910,err=9920) chemps2_lrestart
  ReadStatus = ' O.K. after reading data after DMRE keyword.'
  call ChkIfKey()
end if

!---  Process MXCA command
if (KeyMXCA) then
  call SetPos(LUInput,'MXCA',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading data after MXCA keyword.'
  read(LUInput,*,end=9910,err=9920) max_canonical
  ReadStatus = ' O.K. after reading data after MXCA keyword.'
  call ChkIfKey()
end if

#endif
#endif

!---  Process HFOC command
! This keyword is to specify a user customized orbs occupancies guess.
! It is used by Block and CheMPS2... but it could be useful for other codes.
! Therefore it is now outside the ifdef Block or CheMPS2.
if (KeyHFOC) then
  write(6,*) ' HFOC keyword was given.'
  call SetPos(LUInput,'HFOC',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading after HFOC keyword.'
  write(6,*) 'NASHT, mxact = ',NASHT,mxact
  read(LUInput,*,end=9910,err=9920) (hfocc(i),i=1,NASHT)
  ReadStatus = ' O.K. reading after HFOC keyword.'
  write(6,*) 'HFOCC read in proc_inp of size:',NASHT
  write(6,*) (hfocc(i),i=1,NASHT)
end if

#ifdef _ENABLE_DICE_SHCI_
!---  Process DICE command
if (KeyDICE) then
  DoBlockDMRG = .true.
  write(6,*) 'DICE> (semistochastic) heat bath configuration interaction (SHCI)'
  call SetPos(LUInput,'DICE',Line,iRc)
  call ChkIfKey()
end if
!---  Process STOC command
if (KeySTOC) then
  Dice_Stoc = .true.
  write(6,*) 'DICE> Using semistochastic algorithm interaction (SHCI)'
  call SetPos(LUInput,'STOC',Line,iRc)
  call ChkIfKey()
end if
!---  Process DIOC command
DICEOCC = ''
if (KeyDIOC) then
  call SetPos(LUInput,'DIOC',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading data after DIOC keyword.'
  read(LUInput,*,end=9910,err=9920) nref_dice
  do iref_dice=1,nref_dice
    read(LUInput,'(A)',end=9910,err=9920) diceocc(iref_dice)
    call molcas2dice(diceocc(iref_dice))
  end do
  ReadStatus = ' O.K. after reading data after DIOC keyword.'
  call ChkIfKey()
end if
!---  Process EPSI command
if (KeyEPSI) then
  if (DBG) write(6,*) ' EPS (Thresholds) command was used.'
  call SetPos(LUInput,'EPS',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading thresholds after EPSI keyword.'
  read(LUInput,*,end=9910,err=9920) dice_eps1,dice_eps2
  ReadStatus = ' O.K. after reading thresholds after EPSI keyword.'
  call ChkIfKey()
end if
!---  Process SAMP command
if (KeySAMP) then
  call SetPos(LUInput,'SAMP',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading data after SAMP keyword.'
  read(LUInput,*,end=9910,err=9920) dice_sampleN
  ReadStatus = ' O.K. after reading data after SAMP keyword.'
  call ChkIfKey()
end if
!---  Process DITE command
if (KeyDITE) then
  call SetPos(LUInput,'DITE',Line,iRc)
  if (iRc /= _RC_ALL_IS_WELL_) goto 9810
  ReadStatus = ' Failure reading data after DITE keyword.'
  read(LUInput,*,end=9910,err=9920) dice_iter
  ReadStatus = ' O.K. after reading data after DITE keyword.'
  call ChkIfKey()
end if
!---  Process DIRE command
if (KeyDIRE) then
  dice_restart = .true.
  call SetPos(LUInput,'DIRE',Line,iRc)
  call ChkIfKey()
end if
#endif

!---  All keywords have been processed --------------------------------*

if ((.not. KeyINAC) .and. (sum(nIsh(:nSym)) == 0)) &
  call WarningMessage(1,'The number of inactive orbitals is zero. Do you really want this?')

if (IfCRPR) then
  ! Core shift using a fixed projection operator.
  NCRVEC = NBAS(1)
  call mma_allocate(CRVEC,NCRVEC,Label='CRVec')
  N = NBAS(1)
  NCRPROJ = (N*(N+1)*(N**2+N+2))/8
  call mma_allocate(CRPROJ,NCRPROJ,Label='CRPROJ')
end if
!***********************************************************************
! Generate artificial splitting or RAS into GAS for parallel blocking  *
!***********************************************************************
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
!                                                                      *
!***********************************************************************
!                                                                      *
! Select default root for geometry optimization

if ((NROOTS > 1) .and. (irlxroot == 0)) then

  ! Check if multi state SA-CASSCF

  nW = 0
  do iR=1,LROOTS
    if (WEIGHT(iR) /= 0.0d0) nW = nW+1
  end do
  if (nW /= 1) then
    iRlxRoot = iroot(LROOTS)
  else
    do iR=1,LROOTS
      if (WEIGHT(iR) /= 0.0d0) iRlxRoot = iroot(iR)
    end do
  end if
end if
if ((NROOTS == 1) .or. (LROOTS == 1)) iRlxRoot = iRoot(1)
!                                                                      *
!***********************************************************************
!                                                                      *
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
      do NU=1,NT-1
        ITU = ITU+1
        IZROT(ITU) = 1
      end do
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
  if (irc /= 0) Go To 9930
end if

! ======================================================================

100 continue
!PAM Jump here in case of CASVB (IFVB == 2)
if (DBG) write(6,*) ' After 100 CONTINUE.'

!PAM July 2007 Check in case of CI restart:
if (DBG) write(6,*) ' Check if CI-Restart.'
if (KeyCIRE) then
  ! Test read:
  if (DBG) write(6,*) ' Yes it is!'
  iJOB = -1
  call f_Inquire('JOBOLD',lExists)
  if (lExists) then
    if (DBG) write(6,*) ' ''JOBOLD'' exists.'
    iJOB = 1
  else
    call f_Inquire(IPHNAME,lExists)
    if (lExists) then
      iJOB = 0
      if (DBG) write(6,*) ' No ''JOBOLD'', but JOBIPH exists.'
      if (DBG) write(6,*) ' It is named ',IPHNAME
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
    write(6,*) ' Test read shows that there is no usable interface'
    write(6,*) ' file, necessary for the requested CI restart.'
    write(6,*) ' Most probable reason: the user has forgotten to'
    write(6,*) ' provide this file. The program will continue,'
    write(6,*) ' but there can be no CI restart.'
    ICIRST = 0
  end if
end if
!PAM July 2007 End of addition
! ======================================================================

! Initialize seward

if (DBG) write(6,*) ' Initialize seward.'
nDiff = 0
if (DSCF .or. RF_On() .or. Langevin_On() .or. PCM_On() .or. Do_OFEmb .or. (KSDFT /= 'SCF')) &
  call IniSew(DSCF .or. Langevin_On() .or. PCM_On(),nDiff)
! ======================================================================
#ifdef _DMRG_
domcpdftDMRG = l_casdft .and. doDMRG
twordm_qcm = domcpdftDMRG .or. (.not. KeyCION)
#endif

! Setup part for DMRG calculations
#ifdef _DMRG_
if (keyDMRG .or. doDMRG) then
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

if (DBG) then
  write(6,*) ' Call ChkInp.'
  call XFlush(6)
end if
call ChkInp()
! ======================================================================

! In DMRG-CASSCF, skip GUGA and LUCIA settings
NCONF = 1
if (DoBlockDMRG) goto 9000
! ======================================================================

! Construct the Guga tables

if (.not. (DoNECI .or. Do_CC_CI .or. DumpOnly)) then
  ! right now skip most part of gugactl for GAS, but only call mknsm.
  if (.not. iDoGas) then
    ! DMRG calculation no need the GugaCtl subroutine
#   ifdef _DMRG_
    if (KeyDMRG .or. doDMRG) then
      call mma_deallocate(initial_occ)
      goto 9000
    else
#   endif
      call Timing(Eterna_1,dum1,dum2,dum3)
      if (DBG) write(6,*) ' Call GugaCtl'
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

! Construct the determinant tables

if (DBG) write(6,*) ' Construct the determinant tables.'
MS2 = iSpin-1

! Set variables needed in Lucia_Ini

call iCopy(mxGAS*mxSym,ngssh,1,ngssh_Molcas,1)
call iCopy(mxGAS*2,igsoccx,1,igsoccx_Molcas,1)
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

if (.not. (KeyDMRG .or. DoNECI .or. Do_CC_CI .or. DumpOnly)) then
  ! switch on/off determinants
# ifdef _DMRG_
  if (.not. doDMRG) then
# endif
    ! Initialize LUCIA and determinant control
    call StatusLine('RASSCF: ','Initializing Lucia...')
    call Lucia_Util('Ini')
    ! to get number of CSFs for GAS
    ! and number of determinants to store
    nconf = 0
    nDet = 0
    do i=1,mxsym
      nconf = nconf+ncsasm(i)
      nDet = nDet+ndtasm(i)
    end do
# ifdef _DMRG_
  end if
# endif
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

! ======================================================================
if (ICICH == 1) then
  call mma_allocate(UG2SG_X,NCONF,Label='UG2SG_X')
  call UG2SG(NROOTS,NCONF,NAC,NACTEL,STSYM,IPR,CONF,CFTP,UG2SG_X,ICI,JCJ,CCI,MXROOT)
  call mma_deallocate(UG2SG_X)
end if
! ======================================================================

! faroald initializations
if (DOFARO) then
  if (NSYM > 1) then
    write(6,'(1X,A)') 'FARO keyword was used, but NSYM > 1,'
    write(6,'(1X,A)') 'switching to LUCIA as the CI backend.'
    DOFARO = .false.
  else
    write(6,'(1X,A)') '**EXPERIMENTAL**'
    write(6,'(1X,A)') 'CI backend is FAROALD instead of LUCIA.'
    write(6,'(1X,A)') '**EXPERIMENTAL**'
    call FAROALD_INIT(NACTEL,NASH(1),ISPIN)
    call CITRANS_INIT(NACTEL,NASH(1),ISPIN)
  end if
end if

Go to 9000

!---  Error exits -----------------------------------------------------*
9810 continue
if (IPRLEV >= TERSE) then
  call WarningMessage(2,'Error in input preprocessing.')
  write(6,*) ' PROC_INP: A keyword was found during prescanning'
  write(6,*) ' the input file, but when later trying to locate'
  write(6,*) ' this input, it could not be found. Something has'
  write(6,*) ' happened to the input file, or else there is some'
  write(6,*) ' strange program error.'
  iRc = _RC_INPUT_ERROR_
end if
Go to 9900

9910 continue
call WarningMessage(2,'End of input file during preprocessing.')
call WarningMessage(2,ReadStatus)
if (IPRLEV >= TERSE) write(6,*) ' Error exit 9910 from PROC_INP.'
iRc = _RC_INPUT_ERROR_
Go to 9900

9920 continue
call WarningMessage(2,'Read error during input preprocessing.')
call WarningMessage(2,ReadStatus)
if (IPRLEV >= TERSE) write(6,*) ' Error exit 9920 from PROC_INP.'
iRc = _RC_INPUT_ERROR_
Go to 9900

9930 continue
call WarningMessage(2,'Error during input preprocessing.')
call WarningMessage(2,ReadStatus)
if (IPRLEV >= TERSE) then
  write(6,*) ' Error exit 9930 from PROC_INP.'
  write(6,*) ' Check previous messages in the output'
  write(6,*) ' to find the reason.'
end if
iRc = _RC_INPUT_ERROR_
Go to 9900

!---  Normal exit -----------------------------------------------------*
9000 continue
if (DBG) write(6,*) ' Normal exit from PROC_INP.'

return

!---  Abnormal exit ---------------------------------------------------*
9900 continue
if (DBG) write(6,*) ' Abnormal exit from PROC_INP.'

end subroutine proc_inp
