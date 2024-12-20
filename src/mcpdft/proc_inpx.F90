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
Subroutine Proc_InpX(DSCF,iRc)
  use definitions,only:iwp,wp,u6
  use Fock_util_global,only:DoCholesky
  use Cholesky,only:ChFracMem
  use mcpdft_input,only:mcpdft_options
  use printlevel,only:terse,debug,insane
  use mcpdft_output,only:iPrLoc
  use rasscf_global,only:IPT2,iRoot,lRoots,NAC,NACPAR,NACPR2,NFR,NIN,NO2M,NORBT,NROOTS,NSEC, &
                          nTot3,nTot4,Weight
  use general_data,only:norb,nash,nssh,ndel,nish,nfro,nrs3,nrs2,nrs1,nbas,nconf,nelec3,nhole1, &
                         nactel,stsym,ispin,ntotsp,ntot2,ntot1,nsym,ndelt,invec,jobiph, &
                         jobold,nrs1t,nrs2t,nrs3t,ntot
#ifdef _HDF5_
  Use mh5,Only:mh5_open_file_r,mh5_exists_attr,mh5_exists_dset,mh5_fetch_attr,mh5_fetch_dset,mh5_close_file
  use stdalloc,only:mma_allocate,mma_deallocate
#endif
  implicit none

#include "rasdim.fh"
#include "warnings.h"

  logical,intent(in) :: dscf
  integer(kind=iwp),intent(out) :: irc

  integer(kind=iwp) :: iprlev,IADR19(15),i,iad19
  integer(kind=iwp) :: iorbdata,isym,ndiff
  real(kind=wp) potnucdummy
  logical(kind=iwp) :: DBG,lExists,RunFile_Exists,keyjobi

#ifdef _HDF5_
  character(len=1),allocatable :: typestring(:)
  integer(kind=iwp) :: NBAS_L(8),mh5id,nsym_l
#endif

  integer(kind=iwp),external :: isFreeUnit
  logical(kind=iwp),external :: Langevin_On,PCM_On

  character(len=72) :: ReadStatus,JobTit(mxTit)
  character(len=LENIN8*mxOrb) :: lJobH1
  character(len=2*72) :: lJobH2

  Call StatusLine('MCPDFT:','Processing Input')

!> default for MC-PDFT: read/write from/to JOBIPH-type files
  keyJOBI = .true.
  iRc = _RC_ALL_IS_WELL_

  !> Local print level in this routine:
  IPRLEV = IPRLOC(1)

  DBG = (IPRLEV >= DEBUG)

! ==== Check if there is any runfile ====
  Call F_Inquire('RUNFILE',RunFile_Exists)
  If(DBG) Write(u6,*) ' Inquire about RUNFILE.'
  IF(RunFile_Exists) Then
    If(DBG) Write(u6,*) ' Yes, there is one.'
    NSYM = 0
    Call qpg_iScalar('nSym',lExists)
    IF(lExists) Then
      Call Get_iScalar('nSym',nSym)
      Call Get_iArray('nBas',nBas,nSym)
      If(DBG) Then
        write(u6,*) ' The following information exists on runfile:'
        write(u6,*) ' Nr of symmetries, NSYM:',NSYM
        write(u6,*) ' Nr of basis functions/symmetry:'
        write(u6,'(1x,8I5)')(NBAS(I),I=1,NSYM)
        Call XFlush(6)
      EndIf
    ELSE
      Call WarningMessage(2,'No symmetry info on runfile.')
      write(u6,*) ' There seems to be no information about symmetry'
      write(u6,*) ' on the runfile! This is an unexpected error.'
      Call Quit(_RC_IO_ERROR_READ_)
    ENDIF
  ELSE
    Call WarningMessage(2,'Cannot find runfile.')
    write(u6,*) ' PROC_INP: Cannot find RUNFILE. This is an unexpected error.'
    Call Quit(_RC_IO_ERROR_READ_)
  ENDIF
! ==== End check if there is any runfile ====

! Make these enumerations??
! Also, allow FILE to specify either a binary (JobIph or HDF5 reference for
! example). No reason why we cannot do this.
  iOrbData = 0
! iOrbData=0: no orbital space data is specified
!         >0: specifications from some orbital file (JOBOLD, JOBIPH, HDF5)
  INVEC = 0
! INVEC=0, no source for orbitals (yet)
!       3, take from JOBOLD, or JOBIPH file
!       4, take from an HDF5 file
!       5, take from startorb (instead of jobold or jobiph) NOT IMPLEMENTED

! ---  ==== FILE(ORB) keyword =====
  If(len_trim(mcpdft_options%wfn_file) /= 0) Then
    keyJOBI = .false.
    if(mcpdft_options%is_hdf5_wfn) then
      invec = 4
    else
      invec = 5
      write(u6,*) 'WARNING: cannot specify non-hdf5'
      write(u6,*) 'file with FILE keyword.'
      call abend()
    EndIf
  endif
!---  ==== JOBI(PH) keyword =====
! The following is run, EXCEPT if FILE key is provided an points to an
! HDF5 input file
! I have a feeling that this should only run IF FILE(ORB) key is not passed
  if(keyJOBI) then
    mcpdft_options%wfn_file = "JOBOLD"
    call f_Inquire("JOBOLD",lExists)
    if(.not. lexists) then
      mcpdft_options%wfn_file = "JOBIPH"
      call f_inquire(mcpdft_options%wfn_file,lexists)
      if(.not. lexists) then
        Write(u6,*)
        Write(u6,*) '******************************************'
        Write(u6,*) 'JOBIPH and JOBOLD does not seem to exist, '
        Write(u6,*) 'so the calculation cannot continue.       '
        Write(u6,*) '******************************************'
        Call Abend()
      endif
    endif
    invec = 3

    if(JOBIPH > 0) Then
      Call DaClos(JOBIPH)
      JOBIPH = -1
    endif
    JOBIPH = IsFreeUnit(15)
    CALL DANAME(JOBIPH,mcpdft_options%wfn_file)
    INVEC = 3
  endif !> JOBI(PH) keyword
!---  ==== JOBI(PH) keyword =====

!--- Finish process..some cleanup
  If(mcpdft_options%otfnal%is_hybrid()) Then
    CALL Put_DScalar('R_WF_HMC',mcpdft_options%otfnal%lambda)
    If(DBG) then
      Write(u6,*) 'Wave Funtion Ratio in hybrid PDFT',mcpdft_options%otfnal%lambda
    endif
  EndIf

!---  Process HDF5 file --------------------------------------------*
  If(mcpdft_options%is_hdf5_wfn) Then
#ifdef _HDF5_
    mh5id = mh5_open_file_r(mcpdft_options%wfn_file)
!     read basic attributes
    call mh5_fetch_attr(mh5id,'NSYM',NSYM_L)
    if(nsym /= nsym_l) then
      write(u6,*) 'Number of symmetries on HDF5 file does not'
      write(u6,*) 'match the number of symmetries on the'
      write(u6,*) 'RunFile, calculation will stop now.'
      call Quit(_RC_INPUT_ERROR_)
    endif
    call mh5_fetch_attr(mh5id,'NBAS',NBAS_L)

    do isym = 1,nsym
      if(nbas(isym) /= nbas_l(isym)) then
        write(u6,*) 'Number of basis functions on HDF5 file does not'
        write(u6,*) 'match the number of basis functions on the'
        write(u6,*) 'RunFile, calculation will stop now.'
        call Quit(_RC_INPUT_ERROR_)
      endif
    enddo
!     orbitals available?
    if(.not. mh5_exists_dset(mh5id,'MO_VECTORS')) then
      write(u6,*) 'The HDF5 ref file does not contain MO vectors.'
      write(u6,*) 'Fatal error, the calculation will stop now.'
      call Quit(_RC_INPUT_ERROR_)
    endif
!     typeindex data available?
    if(mh5_exists_dset(mh5id,'MO_TYPEINDICES')) then
      iOrbData = 3
      call mma_allocate(typestring,sum(nbas(1:nsym)))
      call mh5_fetch_dset(mh5id,'MO_TYPEINDICES',typestring)
      call tpstr2orb(nSym,nbas_l,typestring,nfro,nish,NRS1,NRS2,NRS3,nSSh,nDel)
      call mma_deallocate(typestring)
    else
      write(u6,*) 'The HDF5 ref file does not contain TYPEindices.'
      write(u6,*) 'Fatal error, the calculation will stop now.'
      call Quit(_RC_INPUT_ERROR_)
    endif

#ifdef _DMRG_
    if(.not. mh5_exists_dset(mh5id,'QCMAQUIS_CHECKPOINT')) then
#endif
      if(mh5_exists_dset(mh5id,'CI_VECTORS')) then
        write(u6,*) ' CI vectors will be read from HDF5 ref file.'
      else
        write(u6,*) 'The HDF5 ref file does not contain CI vectors.'
        write(u6,*) 'Fatal error, the calculation will stop now.'
        call Quit(_RC_INPUT_ERROR_)
      endif
#ifdef _DMRG_
    endif
#endif

    call mh5_close_file(mh5id)
#else
    write(u6,*) 'The format of the start orbital file was'
    write(u6,*) 'specified by the user as HDF5, but this'
    write(u6,*) 'is not implemented in this installation.'
    call Quit(_RC_INPUT_ERROR_)
#endif
  EndIf

  iprlev = insane

!> read orbital space data AND CI optimiation parameters from JOBIPH
  IF(IORBDATA == 0) THEN
    IAD19 = 0
    Call IDaFile(JOBIPH,2,IADR19,10,IAD19)
    iAd19 = iAdr19(1)
    CALL WR_RASSCF_Info(JobIPH,2,iAd19,NACTEL,ISPIN,NSYM,STSYM,NFRO,NISH,NASH,NDEL,NBAS, &
                        mxSym,lJobH1,LENIN8*mxOrb,NCONF,lJobH2,2*72,JobTit,4*18*mxTit, &
                        POTNUCDUMMY,LROOTS,NROOTS,IROOT,mxRoot,NRS1,NRS2,NRS3,NHOLE1,NELEC3,IPT2,WEIGHT)
  EndIf  !> IORBDATA

!> read CI optimiation parameters from HDF5 file
  if(mcpdft_options%is_hdf5_wfn) then
#ifdef _HDF5_
    mh5id = mh5_open_file_r(mcpdft_options%wfn_file)
    call mh5_fetch_attr(mh5id,'SPINMULT',iSpin)
    call mh5_fetch_attr(mh5id,'NSYM',nSym)
    call mh5_fetch_attr(mh5id,'LSYM',stSym)
    call mh5_fetch_attr(mh5id,'NBAS',nBas)

    call mh5_fetch_attr(mh5id,'NACTEL',nactel)
    call mh5_fetch_attr(mh5id,'NHOLE1',nhole1)
    call mh5_fetch_attr(mh5id,'NELEC3',nelec3)
    call mh5_fetch_attr(mh5id,'NCONF',nconf)
    call mh5_fetch_attr(mh5id,'NSTATES',lroots)
    If(mh5_exists_attr(mh5id,'NROOTS')) Then
      call mh5_fetch_attr(mh5id,'NROOTS',nroots)
    Else
      nroots = lroots
    EndIf
    call mh5_fetch_attr(mh5id,'STATE_WEIGHT',weight)

    call mh5_close_file(mh5id)
#endif
  endif

!AMS - this may be closing either JOBOLD or JOBIPH. Close only JOBOLD.
  IF(JOBOLD > 0) then
    IF(JOBOLD /= JOBIPH) THEN
      Call DaClos(JOBOLD)
    ENDIF
  endif

!AMS - make sure we change to a different JOBIPH file - we don't want to
!overwrite any existing JOBIPH file.
!
!Close the old JOBIPH file
  if(JOBIPH > 0) Then
    Call DaClos(JOBIPH)
    JOBIPH = -1
  endif
!Rename JOBIPH file, and open it.
  JOBIPH = IsFreeUnit(15)
  CALL DANAME(JOBIPH,"JOBIPH")

!---  complete orbital specifications ---------------------------------*
  Do iSym = 1,nSym
    nash(isym) = nrs1(isym)+nrs2(isym)+nrs3(isym)
    NORB(ISYM) = NBAS(ISYM)-NFRO(ISYM)-NDEL(ISYM)
    NSSH(ISYM) = NORB(ISYM)-NISH(ISYM)-NASH(ISYM)
  EndDo
  write(6,*) 'MRH nash:',nash
!---  Related data for sizes, etc.
  NTOT = 0
  NTOT1 = 0
  NTOT2 = 0
  NO2M = 0
  nin = 0
  nac = 0
  NDELT = 0
  nfr = 0
  NSEC = 0
  NORBT = 0
  NTOT3 = 0
  NTOTSP = 0
  NTOT4 = 0
  NRS1T = 0 ! for RASSCF
  NRS2T = 0
  NRS3T = 0
  DO ISYM = 1,NSYM
    NTOT = NTOT+NBAS(ISYM)
    NTOT1 = NTOT1+NBAS(ISYM)*(NBAS(ISYM)+1)/2
    NTOT2 = NTOT2+NBAS(ISYM)**2
    NO2M = MAX(NO2M,NBAS(ISYM)**2)
    NRS1T = NRS1T+NRS1(ISYM)  ! for RAS
    NRS2T = NRS2T+NRS2(ISYM)
    NRS3T = NRS3T+NRS3(ISYM)
    nfr = nfr+NFRO(ISYM)
    nin = nin+NISH(ISYM)
    nac = nac+NASH(ISYM)
    NDELT = NDELT+NDEL(ISYM)
    NSEC = NSEC+NSSH(ISYM)
    NORBT = NORBT+NORB(ISYM)
    NTOT3 = NTOT3+(NORB(ISYM)+NORB(ISYM)**2)/2
    NTOTSP = NTOTSP+(NASH(ISYM)*(NASH(ISYM)+1)/2)
    NTOT4 = NTOT4+NORB(ISYM)**2
  ENDDO
  NACPAR = (nac+nac**2)/2
  NACPR2 = (NACPAR+NACPAR**2)/2
  write(6,*) 'MRH nac:',nac

  ! Initialize Cholesky information if requested
  if(DoCholesky) then
    Call Cho_X_init(irc,ChFracMem)
    if(irc /= 0) Go To 9930
  endif

! ===============================================================
!     Initialize seward
  If(DBG) write(u6,*) ' Initialize seward.'
  nDiff = 0
  Call IniSew(DSCF .or. Langevin_On() .or. PCM_On(),nDiff)
! ===============================================================
!     Check the input data
  write(6,*) 'seward nac:',nac
  Call validate_wfn()
! ===============================================================
  Go to 9000

!---  Error exits -----------------------------------------------------*
9930 CONTINUE
  Call WarningMessage(2,'Error during input preprocessing.')
  Call WarningMessage(2,ReadStatus)
  If(IPRLEV >= TERSE) write(u6,*) ' Error exit 9930 from PROC_INP.'
  iRc = _RC_INPUT_ERROR_
  Go to 9900

!---  Normal exit -----------------------------------------------------*
9000 CONTINUE
  close(989)
  If(DBG) write(u6,*) ' Normal exit from PROC_INP.'
  Return
!---  Abnormal exit -----------------------------------------------------*
9900 CONTINUE
  If(DBG) write(u6,*) ' Abnormal exit from PROC_INP.'
  Return
EndSubroutine
