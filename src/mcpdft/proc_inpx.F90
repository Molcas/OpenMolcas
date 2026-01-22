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

subroutine Proc_InpX(DSCF,iRc)

use Index_Functions, only: nTri_Elem
use Fock_util_global, only: DoCholesky
use Cholesky, only: ChFracMem
use PrintLevel, only: DEBUG, INSANE, TERSE
use mcpdft_input, only: mcpdft_options
use mcpdft_output, only: iPrLoc
use rasscf_global, only: IPT2, iRoot, lRoots, NAC, NACPAR, NACPR2, NFR, NIN, NO2M, NORBT, NROOTS, NSEC, nTot3, nTot4, Weight
use general_data, only: invec, ispin, jobiph, jobold, nactel, nash, nbas, nconf, ndel, ndelt, nelec3, nfro, nhole1, nish, norb, &
                        nrs1, nrs1t, nrs2, nrs2t, nrs3, nrs3t, nssh, nsym, ntot, ntot1, ntot2, ntotsp, stsym
use rctfld_module, only: lrf
#ifdef _HDF5_
use mh5, only: mh5_close_file, mh5_exists_attr, mh5_exists_dset, mh5_fetch_attr, mh5_fetch_dset, mh5_open_file_r
use stdalloc, only: mma_allocate, mma_deallocate
#endif
use Molcas, only: LenIn8, MxOrb, MxRoot, MxSym
use RASDim, only: MxTit
use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp), intent(in) :: dscf
integer(kind=iwp), intent(out) :: irc
#include "warnings.h"
integer(kind=iwp) :: i, iad19, IADR19(15), iorbdata, iprlev, isym, ndiff
real(kind=wp) :: potnucdummy
logical(kind=iwp) :: DBG, keyjobi, lExists, RunFile_Exists
character(len=LenIn8*mxOrb) :: lJobH1
character(len=72) :: JobTit(mxTit), ReadStatus
character(len=2*72) :: lJobH2
#ifdef _HDF5_
integer(kind=iwp) :: mh5id, NBAS_L(MxSym), nsym_l
character, allocatable :: typestring(:)
#endif
integer(kind=iwp), external :: isFreeUnit
logical(kind=iwp), external :: Langevin_On, PCM_On

call StatusLine('MCPDFT:','Processing Input')

!> default for MC-PDFT: read/write from/to JOBIPH-type files
keyJOBI = .true.
iRc = _RC_ALL_IS_WELL_

!> Local print level in this routine:
IPRLEV = IPRLOC(1)

DBG = (IPRLEV >= DEBUG)

! ==== Check if there is any runfile ====
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
if (len_trim(mcpdft_options%wfn_file) /= 0) then
  keyJOBI = .false.
  if (mcpdft_options%is_hdf5_wfn) then
    invec = 4
  else
    invec = 5
    write(u6,*) 'WARNING: cannot specify non-hdf5'
    write(u6,*) 'file with FILE keyword.'
    call abend()
  end if
end if
!---  ==== JOBI(PH) keyword =====
! The following is run, EXCEPT if FILE key is provided an points to an
! HDF5 input file
! I have a feeling that this should only run IF FILE(ORB) key is not passed
if (keyJOBI) then
  mcpdft_options%wfn_file = 'JOBOLD'
  call f_Inquire('JOBOLD',lExists)
  if (.not. lexists) then
    mcpdft_options%wfn_file = 'JOBIPH'
    call f_inquire(mcpdft_options%wfn_file,lexists)
    if (.not. lexists) then
      write(u6,*)
      write(u6,*) '******************************************'
      write(u6,*) 'JOBIPH and JOBOLD does not seem to exist, '
      write(u6,*) 'so the calculation cannot continue.       '
      write(u6,*) '******************************************'
      call Abend()
    end if
  end if
  invec = 3

  if (JOBIPH > 0) then
    call DaClos(JOBIPH)
    JOBIPH = -1
  end if
  JOBIPH = IsFreeUnit(15)
  call DANAME(JOBIPH,mcpdft_options%wfn_file)
  INVEC = 3
end if !> JOBI(PH) keyword
!---  ==== JOBI(PH) keyword =====

!--- Finish process..some cleanup
if (mcpdft_options%otfnal%is_hybrid()) then
  call Put_DScalar('R_WF_HMC',mcpdft_options%otfnal%lambda)
  if (DBG) write(u6,*) 'Wave Funtion Ratio in hybrid PDFT',mcpdft_options%otfnal%lambda
end if

!---  Process HDF5 file -----------------------------------------------*
if (mcpdft_options%is_hdf5_wfn) then
# ifdef _HDF5_
  mh5id = mh5_open_file_r(mcpdft_options%wfn_file)
  ! read basic attributes
  call mh5_fetch_attr(mh5id,'NSYM',NSYM_L)
  if (nsym /= nsym_l) then
    write(u6,*) 'Number of symmetries on HDF5 file does not'
    write(u6,*) 'match the number of symmetries on the'
    write(u6,*) 'RunFile, calculation will stop now.'
    call Quit(_RC_INPUT_ERROR_)
  end if
  call mh5_fetch_attr(mh5id,'NBAS',NBAS_L)

  do isym=1,nsym
    if (nbas(isym) /= nbas_l(isym)) then
      write(u6,*) 'Number of basis functions on HDF5 file does not'
      write(u6,*) 'match the number of basis functions on the'
      write(u6,*) 'RunFile, calculation will stop now.'
      call Quit(_RC_INPUT_ERROR_)
    end if
  end do
  ! orbitals available?
  if (.not. mh5_exists_dset(mh5id,'MO_VECTORS')) then
    write(u6,*) 'The HDF5 ref file does not contain MO vectors.'
    write(u6,*) 'Fatal error, the calculation will stop now.'
    call Quit(_RC_INPUT_ERROR_)
  end if
  ! typeindex data available?
  if (mh5_exists_dset(mh5id,'MO_TYPEINDICES')) then
    iOrbData = 3
    call mma_allocate(typestring,sum(nbas(1:nsym)))
    call mh5_fetch_dset(mh5id,'MO_TYPEINDICES',typestring)
    call tpstr2orb(nSym,nbas_l,typestring,nfro,nish,NRS1,NRS2,NRS3,nSSh,nDel)
    call mma_deallocate(typestring)
  else
    write(u6,*) 'The HDF5 ref file does not contain TYPEindices.'
    write(u6,*) 'Fatal error, the calculation will stop now.'
    call Quit(_RC_INPUT_ERROR_)
  end if

# ifdef _DMRG_
  if (.not. mh5_exists_dset(mh5id,'QCMAQUIS_CHECKPOINT')) then
# endif
    if (mh5_exists_dset(mh5id,'CI_VECTORS')) then
      write(u6,*) ' CI vectors will be read from HDF5 ref file.'
    else
      write(u6,*) 'The HDF5 ref file does not contain CI vectors.'
      write(u6,*) 'Fatal error, the calculation will stop now.'
      call Quit(_RC_INPUT_ERROR_)
    end if
# ifdef _DMRG_
  end if
# endif

  call mh5_close_file(mh5id)
# else
  write(u6,*) 'The format of the start orbital file was'
  write(u6,*) 'specified by the user as HDF5, but this'
  write(u6,*) 'is not implemented in this installation.'
  call Quit(_RC_INPUT_ERROR_)
# endif
end if

iprlev = INSANE

!> read orbital space data AND CI optimiation parameters from JOBIPH
if (IORBDATA == 0) then
  IAD19 = 0
  call IDaFile(JOBIPH,2,IADR19,10,IAD19)
  iAd19 = iAdr19(1)
  call WR_RASSCF_Info(JobIPH,2,iAd19,NACTEL,ISPIN,NSYM,STSYM,NFRO,NISH,NASH,NDEL,NBAS,mxSym,lJobH1,LENIN8*mxOrb,NCONF,lJobH2,2*72, &
                      JobTit,4*18*mxTit,POTNUCDUMMY,LROOTS,NROOTS,IROOT,mxRoot,NRS1,NRS2,NRS3,NHOLE1,NELEC3,IPT2,WEIGHT)
end if  !> IORBDATA

!> read CI optimization parameters from HDF5 file
#ifdef _HDF5_
if (mcpdft_options%is_hdf5_wfn) then
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
  if (mh5_exists_attr(mh5id,'NROOTS')) then
    call mh5_fetch_attr(mh5id,'NROOTS',nroots)
  else
    nroots = lroots
  end if
  call mh5_fetch_attr(mh5id,'STATE_WEIGHT',weight)

  call mh5_close_file(mh5id)
end if
#endif

! AMS - this may be closing either JOBOLD or JOBIPH. Close only JOBOLD.
if (JOBOLD > 0) then
  if (JOBOLD /= JOBIPH) call DaClos(JOBOLD)
end if

! AMS - make sure we change to a different JOBIPH file - we don't want to
! overwrite any existing JOBIPH file.
!
!Close the old JOBIPH file
if (JOBIPH > 0) then
  call DaClos(JOBIPH)
  JOBIPH = -1
end if
! Rename JOBIPH file, and open it.
JOBIPH = IsFreeUnit(15)
call DANAME(JOBIPH,'JOBIPH')

!---  complete orbital specifications ---------------------------------*
do iSym=1,nSym
  nash(isym) = nrs1(isym)+nrs2(isym)+nrs3(isym)
  NORB(ISYM) = NBAS(ISYM)-NFRO(ISYM)-NDEL(ISYM)
  NSSH(ISYM) = NORB(ISYM)-NISH(ISYM)-NASH(ISYM)
end do
!---  Related data for sizes, etc.
NTOT1 = 0
NTOT2 = 0
NO2M = 0
NTOT3 = 0
NTOTSP = 0
NTOT4 = 0
ntot = sum(nbas(1:nsym))
nrs1t = sum(nrs1(1:nsym)) ! for RAS
nrs2t = sum(nrs2(1:nsym))
nrs3t = sum(nrs3(1:nsym))
nfr = sum(nfro(1:nsym))
nin = sum(nish(1:nsym))
nac = sum(nash(1:nsym))
ndelt = sum(ndel(1:nsym))
nsec = sum(nssh(1:nsym))
norbt = sum(norb(1:nsym))
do ISYM=1,NSYM
  NTOT1 = NTOT1+nTri_Elem(NBAS(ISYM))
  NTOT2 = NTOT2+NBAS(ISYM)**2
  NO2M = max(NO2M,NBAS(ISYM)**2)
  ntot3 = ntot3+nTri_Elem(norb(isym))
  NTOTSP = NTOTSP+nTri_Elem(NASH(ISYM))
  NTOT4 = NTOT4+NORB(ISYM)**2
end do
NACPAR = nTri_Elem(nac)
NACPR2 = nTri_Elem(NACPAR)

! Initialize Cholesky information if requested
if (DoCholesky) then
  call Cho_X_init(irc,ChFracMem)
  if (irc /= 0) then
    !---  Error exit --------------------------------------------------*
    call WarningMessage(2,'Error during input preprocessing.')
    call WarningMessage(2,ReadStatus)
    if (IPRLEV >= TERSE) write(u6,*) ' Error exit from PROC_INP.'
    iRc = _RC_INPUT_ERROR_
    if (DBG) write(u6,*) ' Abnormal exit from PROC_INP.'
    return
  end if
end if

! ===============================================================
!     Initialize seward
if (DBG) write(u6,*) ' Initialize seward.'
nDiff = 0
call IniSew(DSCF .or. Langevin_On() .or. PCM_On(),nDiff)
if (lrf) then
  call warningmessage(2,'MC-PDFT with solvent not supported!')
  write(u6,*) 'MC-PDFT cannot be used with solvent!'
  call Quit_OnUserError()
end if
! ===============================================================
!     Check the input data
call validate_wfn()
! ===============================================================
!---  Normal exit -----------------------------------------------------*
if (DBG) write(u6,*) ' Normal exit from PROC_INP.'

end subroutine Proc_InpX
