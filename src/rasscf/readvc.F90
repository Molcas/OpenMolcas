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
! Copyright (C) 1998, Markus P. Fuelscher                              *
!               2018, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine ReadVC(CMO,OCC,D,DS,P,PA,scheme)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     get start MO-coefficients                                        *
!     (if the CI is restarted from a previous calculation get          *
!      also the density matrices)                                      *
!                                                                      *
!     calling arguments:                                               *
!     CMO     : array of real*8                                        *
!               MO-coefficients                                        *
!     OCC     : array of real*8                                        *
!               occupation numbers                                     *
!     INVEC   : integer                                                *
!               flag indicating orbital type                           *
!     D       : array of real*8                                        *
!               averaged one-body density matrix                       *
!     DS      : array of real*8                                        *
!               averaged one-body spin density matrix                  *
!     P       : array of real*8                                        *
!               averaged two body density matrix                       *
!     PA      : array of real*8                                        *
!               averaged antisymmetric twobody density matrix          *
!     FI      : array of real*8                                        *
!               inactive Fock matrix                                   *
!     D1I     : array of real*8                                        *
!               inactive one body density matrix                       *
!     D1A     : array of real*8                                        *
!               active one body density matrix                         *
!     TUVX    : array of real*8                                        *
!               two-electron integrals (tu!vx)                         *
!     IFINAL  : integer                                                *
!               termination flag                                       *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M. P. Fuelscher                                                  *
!     University of Lund, Sweden, 1998                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history:                                                         *
!     none                                                             *
!                                                                      *
!***********************************************************************

use rasscf_global, only: lRoots, nRoots, iRoot, Weight, nAcPar, iXsym, iAlphaBeta, iOverwr, iSUPSM, iCIrst, iPhName, nAcpr2, &
                         nOrbT, purify, iAdr15
use general_data, only: nSym, nDel, nBas, nOrb, nTot, nTot2, Invec, LuStartOrb, StartOrbFile, JobOld, JobIph, nSSH
use casvb_global, only: ifvb
use orthonormalization, only: t_ON_scheme, ON_scheme_values, orthonormalize
use general_data, only: CleanMask
use PrintLevel, only: DEBUG, TERSE, VERBOSE
use output_ras, only: LF, IPRGLB, IPRLOC
use Molcas, only: LenIn, MaxBfn, MxOrb, MxRoot, MxSym
use RASDim, only: MxTit
#ifdef _HDF5_
use mh5, only: mh5_open_file_r, mh5_exists_dset, mh5_fetch_dset, mh5_close_file
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: RtoI

implicit none
character(len=16), parameter :: ROUTINE = 'READVC  '
real*8 :: CMO(*), OCC(*), D(*), DS(*), P(*), PA(*)
type(t_ON_scheme), intent(in) :: scheme
logical :: found, changed
integer :: iPrlev, nData, i, j, NNwOrd, iSym, iErr, IAD19, iJOB, lll, iDisk, jRoot, kRoot, iDummy(1), IADR19(30), iAD15, nTmp(8)
real*8 :: Dummy(1), Scal
real*8, allocatable :: CMO_copy(:)
character(len=(LenIn+8)*mxOrb) :: lJobH1
character(len=2*72) :: lJobH2
character(len=72) :: JobTit(mxTit)
character(len=80) :: VecTit
character(len=4) :: Label
integer, allocatable :: TIND(:), NewOrd(:), TmpXSym(:), JobH(:)
real*8, allocatable :: Scr(:), Ene(:), JobR(:)
#ifdef _HDF5_
integer mh5id
character(len=maxbfn) typestring
#endif
interface
  integer function isfreeunit(seed)
    integer, intent(in) :: seed
  end function
end interface
#include "warnings.h"

!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
! Local print level (if any)
IPRLEV = IPRLOC(1)
if (IPRLEV >= DEBUG) write(LF,*) ' Entering ',ROUTINE
!----------------------------------------------------------------------*
! Do we use default orbitals?                                          *
!----------------------------------------------------------------------*
if (InVec == 0) then
  call qpg_darray('RASSCF orbitals',Found,nData)
  if (Found) then
    InVec = 6
    if (IPRLEV >= TERSE) write(6,'(6x,a)') 'Orbitals from runfile: rasscf orbitals'
  end if
end if
call Check_InVec(InVec)
if (InVec == 0) then
  call qpg_darray('SCF orbitals',Found,nData)
  if (Found) then
    InVec = 7
    if (IPRLEV >= TERSE) write(6,'(6x,a)') 'Orbitals from runfile: scf orbitals'
  end if
end if
call Check_InVec(InVec)
if (InVec == 0) then
  call qpg_darray('Guessorb',Found,nData)
  if (Found) then
    InVec = 5
    if (IPRLEV >= TERSE) write(6,'(6x,a)') 'Orbitals from runfile: guessorb orbitals'
  end if
end if
call Check_InVec(InVec)
if (Invec == 0) InVec = 1

!----------------------------------------------------------------------*
LUStartOrb = 19
LUStartOrb = IsFreeUnit(LUStartOrb)
if (ifvb == 2) invec = 3
if (InVec == 2) then
  ! read from unit formatted ascii file with starting orbitals

  ! Note: Inside RDVEC, the file StartOrbFile is opened, but uses blindly
  ! the unit number provided here. So that should better be a usable
  ! number, or else!
  Label = 'CO  '
  if (iAlphaBeta == 1) Label(3:3) = 'A'
  if (iAlphaBeta == -1) Label(3:3) = 'B'
  if (iOverwr == 1) then
    call RDVEC(StartOrbFile,LUStartOrb,Label,NSYM,NBAS,NBAS,CMO,OCC,Dummy,iDummy,VECTIT,0,iErr)
  else
    Label(4:4) = 'I'
    call mma_allocate(TIND,maxbfn,Label='TIND')
    call RDVEC(StartOrbFile,LUStartOrb,Label,NSYM,NBAS,NBAS,CMO,OCC,Dummy,TIND,VECTIT,0,iErr)
    ! If the typeindex array is used to resort orbitals, then if also
    ! a supersymmetry array is used, it has to be changed.
    ! The supersymmtry array is IXSYM().
    ! VECSORT is a utility that does not know about supersymmetry.
    ! So changing any orbital indices in IXSYM (or potentially any
    ! other orbital indices -- what about ALTER??) must be done HERE
    ! immediately togather with the VecSort, but not *inside* VecSort.

    ! But VecSort does not return any indexing information -- how are
    ! we to know how to change IXSYM?
    ! VecSort changed to include a reindexing array!
    NNwOrd = 0
    do ISym=1,NSym
      NNwOrd = NNwOrd+NBas(ISym)
    end do
    call mma_allocate(NEWORD,NNwOrd,Label='NewOrd')
    !call VecSort(NSYM,NBAS,NBAS,CMO,OCC,TIND,iErr)
    call VecSort(NSYM,NBAS,NBAS,CMO,OCC,TIND,NNwOrd,NewOrd,iErr)
    ! If there is a supersymmetry array, use the orbital mapping:
    if (iSUPSM /= 0) then
      call mma_allocate(TMPXSYM,NNwOrd,Label='TmpXSym')
      do I=1,NNwOrd
        J = NewOrd(I)
        TmpXSym(I) = IXSYM(J)
      end do
      call ICopy(NNwOrd,TmpXSym,1,IXSYM,1)
      call mma_deallocate(TMPXSYM)
    end if

    call mma_deallocate(NewOrd)
    call mma_deallocate(TIND)
  end if
  close(LUStartOrb)
  if (iErr == 1) then
    write(LF,*) 'RASSCF tried to read input orbitals from a'
    write(LF,*) 'file, but encountered an error in the'
    write(LF,*) 'TypeIndex data.'
    call Abend()
    return
  end if

  if (IPRLEV >= TERSE) then
    write(LF,'(6X,A)') 'The MO-coefficients are taken from the file:'
    write(LF,'(6X,A)') trim(StartOrbFile)
    write(LF,'(6X,A,A)') 'Title:',VecTit(2:80)
  end if

else if (InVec == 3) then
  ! read from unit JOBOLD (binary file)
  IAD19 = 0
  iJOB = 0
  call f_Inquire('JOBOLD',Found)
  if (Found) iJOB = 1
  if (iJOB == 1) then
    if (JOBOLD <= 0) then
      JOBOLD = 20
      call DaName(JOBOLD,'JOBOLD')
    end if
  else
    if (IPRLEV >= TERSE) write(LF,'(6X,A)') 'File JOBOLD not found -- use JOBIPH.'
    if (JOBIPH > 0) then
      JOBOLD = JOBIPH
    else
      call DaName(JOBOLD,IPHNAME)
    end if
  end if
  call IDaFile(JOBOLD,2,IADR19,15,IAD19)
  if (IADR19(15) == -1) then
    IAD19 = 0
    call IDAFILE(JOBOLD,2,IADR19,30,IAD19)
  else
    do I=16,30
      IADR19(I) = 0
    end do
    if (IPRGLB >= VERBOSE) call WarningMessage(1,'Old JOBIP file layout.')
  end if
  lll = 10+RtoI
  lll = max(lll,mxSym)
  lll = max(lll,mxOrb)
  lll = max(lll,mxRoot)
  call mma_allocate(JobH,lll,Label='JobH')
  call mma_allocate(JobR,MxRoot,Label='JobH')
  iAd19 = iAdr19(1)
  call WR_RASSCF_Info(JobOld,2,iAd19,JobH(1),JobH(2),JobH(3),JobH(4),JobH,JobH,JobH,JobH,JobH,mxSym,lJobH1,(LenIn+8)*mxOrb, &
                      JobH(5),lJobH2,2*72,JobTit,72*mxTit,JobR(1),JobH(6),JobH(7),JobH,mxRoot,JobH,JobH,JobH,JobH(8),JobH(9), &
                      JobH(10),JobR(:))
  if (IPRLEV >= TERSE) then
    if (iJOB == 1) then
      write(LF,'(6X,A)') 'The MO-coefficients are taken from the file:'
      write(LF,'(6X,A)') 'JOBOLD'
    else
      write(LF,'(6X,A)') 'The MO-coefficients are taken from the file:'
      write(LF,'(6X,A)') trim(iPhName)
    end if
    write(VecTit(1:72),'(A72)') JobTit(1)
    write(LF,'(6X,2A)') 'Title:',VecTit(1:72)
  end if
  call mma_deallocate(JobH)
  call mma_deallocate(JobR)
  iAd19 = iAdr19(2)
  call DDaFile(JobOld,2,CMO,NTOT2,iAd19)
  call DDaFile(JobOld,2,OCC,nTot,iAd19)
  if (ICIRST == 1) then
    if (IPRLEV >= VERBOSE) then
      if (iJOB == 1) then
        write(LF,'(6X,A)') 'The active density matrices (D,DS,P,PA) are read from file JOBOLD and weighted together.'
      else
        write(LF,'(6X,A)') 'The active density matrices (D,DS,P,PA) are read from file '//trim(iPhName)//' and weighted together.'
      end if
    end if
    call mma_allocate(Scr,NACPR2,Label='Scr')
    iDisk = IADR19(3)
    do jRoot=1,lRoots
      Scal = 0.0d0
      do kRoot=1,nRoots
        if (iRoot(kRoot) == jRoot) Scal = Weight(kRoot)
      end do
      call DDaFile(JOBOLD,2,scr,NACPAR,iDisk)
      call daxpy_(NACPAR,Scal,scr,1,D,1)
      call DDaFile(JOBOLD,2,scr,NACPAR,iDisk)
      call daxpy_(NACPAR,Scal,scr,1,DS,1)
      call DDaFile(JOBOLD,2,scr,NACPR2,iDisk)
      call daxpy_(NACPR2,Scal,scr,1,P,1)
      call DDaFile(JOBOLD,2,scr,NACPR2,iDisk)
      call daxpy_(NACPR2,Scal,scr,1,PA,1)
    end do
    call mma_deallocate(Scr)
  end if
  !SVC: read the L2ACT and LEVEL arrays from the jobiph file
  !IFG: disabled, since it breaks when changing active space specification
  !IAD19 = IADR19(18)
  !if (IAD19 /= 0) then
  !  call IDAFILE(JOBOLD,2,IDXSX,mxAct,IAD19)
  !  call IDAFILE(JOBOLD,2,IDXCI,mxAct,IAD19)
  !end if
  if ((JOBOLD > 0) .and. (JOBOLD /= JOBIPH)) then
    call DaClos(JOBOLD)
    JOBOLD = -1
  else if (JOBOLD > 0) then
    JOBOLD = -1
  end if

else if (InVec == 4) then
  ! read from a HDF5 wavefunction file
# ifdef _HDF5_
  if (IPRLEV >= TERSE) then
    write(LF,'(6X,A)') 'The MO-coefficients are taken from the file:'
    write(LF,'(6X,A)') trim(StartOrbFile)
  end if

  mh5id = mh5_open_file_r(StartOrbFile)
  typestring = ''
  select case (iAlphaBeta)
    case (1)
      Label = 'CA  '
      VecTit = 'MO_ALPHA_TYPEINDICES'
    case (-1)
      Label = 'CB  '
      VecTit = 'MO_BETA_TYPEINDICES'
    case default
      Label = 'C   '
      VecTit = 'MO_TYPEINDICES'
  end select
  call RdVec_HDF5(mh5id,Label,NSYM,NBAS,CMO,Dummy,Dummy,iDummy)
  if (mh5_exists_dset(mh5id,trim(VecTit))) call mh5_fetch_dset(mh5id,trim(VecTit),typestring)
  call mh5_close_file(mh5id)
  ! Reorder orbitals based on typeindex
  if (typestring /= '') then
    NNwOrd = 0
    do iSym=1,nSym
      NNwOrd = NNwOrd+NBas(iSym)
    end do
    call mma_allocate(TInd,NNWOrd,Label='TIND')
    call mma_allocate(NewOrd,NNwOrd,Label='NewOrd')
    call tpstr2tpidx(typestring,TInd,NNWOrd)
    call VecSort(NSYM,NBAS,NBAS,CMO,OCC,TInd,NNwOrd,NewOrd,iErr)
    ! If there is a supersymmetry array, use the orbital mapping:
    if (iSUPSM /= 0) then
      call mma_allocate(TmpXSym,NNwOrd,Label='TmpXSym')
      do i=1,NNwOrd
        j = NewOrd(i)
        TmpXSym(i) = iXSym(j)
      end do
      call iCopy(NNwOrd,TmpXSym,1,iXSym,1)
      call mma_deallocate(TmpXSym)
    end if
    call mma_deallocate(NewOrd)
    call mma_deallocate(TInd)
  end if
# else
  write(6,*) 'Orbitals requested from HDF5, but this'
  write(6,*) 'installation does not support that, abort!'
  call abend()
# endif

  ! guess MO-coefficients

else if (InVec == 5) then
  if (IPRLEV >= VERBOSE) write(LF,'(6x,a)') 'Detected guessorb orbitals'
  call Qpg_dArray('Guessorb',Found,nData)
  call Get_dArray('Guessorb',CMO,nData)
  call Qpg_iArray('nDel_go',Found,nData)
  if (Found) then
    call Get_iArray('nDel_go',nTmp,nData)
    Changed = .false.
    do iSym=1,nSym
      if (nTmp(iSym) > nDel(iSym)) Changed = .true.
    end do
    if (Changed) then
      write(6,'(5x,a,8i5)') 'Number of deleted orbitals changed from',(nDel(i),i=1,nSym)
      write(6,'(5x,a,8i5)') '                           changed to  ',(nTmp(i),i=1,nSym)
    end if
    do iSym=1,nSym
      if (nTmp(iSym) > nDel(iSym)) then
        nSsh(iSym) = nSsh(iSym)-nTmp(iSym)+nDel(iSym)
        nOrb(iSym) = nOrb(iSym)-nTmp(iSym)+nDel(iSym)
        nOrbT = nOrbT-nTmp(iSym)+nDel(iSym)
        nDel(iSym) = nTmp(iSym)
      end if
    end do
    if (IPRLEV >= TERSE) write(LF,'(6X,A)') 'The MO-coefficients are taken from guessorb on runfile'
  end if
else if (InVec == 6) then
  if (IPRLEV >= VERBOSE) write(LF,'(6x,a)') 'Detected old RASSCF orbitals'
  call qpg_darray('RASSCF orbitals',Found,nData)
  call get_darray('RASSCF orbitals',CMO,nData)
  call Qpg_iArray('nDel',Found,nData)
  if (Found) then
    call Get_iArray('nDel',nTmp,nData)
    Changed = .false.
    do iSym=1,nSym
      if (nTmp(iSym) > nDel(iSym)) Changed = .true.
    end do
    if (Changed) then
      write(6,'(5x,a,8i5)') 'Number of deleted orbitals changed from',(nDel(i),i=1,nSym)
      write(6,'(5x,a,8i5)') '                           changed to  ',(nTmp(i),i=1,nSym)
    end if
    do iSym=1,nSym
      if (nTmp(iSym) > nDel(iSym)) then
        nSsh(iSym) = nSsh(iSym)-nTmp(iSym)+nDel(iSym)
        nOrb(iSym) = nOrb(iSym)-nTmp(iSym)+nDel(iSym)
        nOrbT = nOrbT-nTmp(iSym)+nDel(iSym)
        nDel(iSym) = nTmp(iSym)
      end if
    end do
    if (IPRLEV >= TERSE) write(LF,'(6X,A)') 'The MO-coefficients are taken from rasscf orbitals on runfile'
  end if
else if (InVec == 7) then
  if (IPRLEV >= VERBOSE) write(LF,'(6x,a)') 'Detected SCF orbitals'
  call qpg_darray('SCF orbitals',Found,nData)
  call get_darray('SCF orbitals',CMO,nData)
  call Qpg_iArray('nDel',Found,nData)
  if (Found) then
    call Get_iArray('nDel',nTmp,nData)
    Changed = .false.
    do iSym=1,nSym
      if (nTmp(iSym) > nDel(iSym)) Changed = .true.
    end do
    if (Changed) then
      write(6,'(5x,a,8i5)') 'Number of deleted orbitals changed from',(nDel(i),i=1,nSym)
      write(6,'(5x,a,8i5)') '                           changed to  ',(nTmp(i),i=1,nSym)
    end if
    do iSym=1,nSym
      if (nTmp(iSym) > nDel(iSym)) then
        nSsh(iSym) = nSsh(iSym)-nTmp(iSym)+nDel(iSym)
        nOrb(iSym) = nOrb(iSym)-nTmp(iSym)+nDel(iSym)
        nOrbT = nOrbT-nTmp(iSym)+nDel(iSym)
        nDel(iSym) = nTmp(iSym)
      end if
    end do
  end if
  if (IPRLEV >= TERSE) write(LF,'(6X,A)') 'The MO-coefficients are taken from scf orbitals on runfile'
else if (InVec == 1) then
  if (IPRLEV >= VERBOSE) write(LF,'(6X,A)') 'The MO-coefficients are obtained by diagonalizing the core Hamiltonian'
  call Guess(CMO)
else
  write(LF,*) 'Severe internal bug prevents further calculation.'
  write(LF,*) 'Invalid value for INVEC in READVC. Program stops.'
  write(LF,*) 'Please issue bug report. INVEC=',INVEC
  call QUIT(_RC_GENERAL_ERROR_)
end if
! print start orbitals
if (IPRLEV >= DEBUG) then
  call mma_allocate(ENE,nTot,Label='ENE')
  call DCOPY_(nTot,[0.0d0],0,ENE,1)
  call PRIMO_RASSCF('Input orbitals',ENE,OCC,CMO)
  call mma_deallocate(ENE)
end if

! cleaning orbitals for high symmetry cases

if (allocated(CleanMask)) call ClnMO(CMO)
if (PURIFY(1:6) == 'LINEAR') call LINPUR(CMO)
if (PURIFY(1:4) == 'ATOM') call SPHPUR(CMO)

if (scheme%val /= ON_scheme_values%no_ON) then
  call mma_allocate(CMO_copy,nTot2)
  CMO_copy(:nTot2) = CMO(:nTot2)
  call orthonormalize(CMO_copy,scheme,CMO(:nTot2))
  call mma_deallocate(CMO_copy)
end if

! save start orbitals

IAD15 = IADR15(2)
call DDAFILE(JOBIPH,1,CMO,NTOT2,IAD15)
call DDAFILE(JOBIPH,1,OCC,NTOT,IAD15)

! exit

return

end subroutine ReadVC
