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

module write_orbital_files

use general_data, only: iSpin, nActel, nAsh, nBas, nConf, nDel, nElec3, nFro, nHole1, nIsh, nRs1, nRs2, nRs3, nSym, nTot, nTot2, &
                        stSym
use gas_data, only: iDoGas, nGAS, nGssh
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
private

logical :: write_orb_per_iter = .false.

interface get_typeidx
  module procedure :: RAS_get_typeidx, GAS_get_typeidx
end interface get_typeidx

interface
  function isfreeunit(iseed)
    import :: iwp
    integer(kind=iwp) :: isfreeunit
    integer(kind=iwp), intent(in) :: iseed
  end function
end interface

public :: get_typeidx, OrbFiles, putOrbFile, write_orb_per_iter

contains

subroutine OrbFiles(JOBIPH,IPRLEV)

  use fortran_strings, only: str
  use rasscf_global, only: BName, E2Act, FDiag, header, iOrbTyp, iPt2, iRoot, iToc, lRoots, maxorbout, nRoots, title, Weight
  use PrintLevel, only: USUAL
  use Molcas, only: LenIn, MxOrb, MxRoot, MxSym
  use RASDim, only: MxIter, MxTit
  use Definitions, only: u6

  implicit none
  integer(kind=iwp), intent(in) :: JobIph, iPrlev
  integer(kind=iwp) :: iDisk, iNDType(7,8), iRt, lUVVVec
  real(kind=wp) :: Energy, PotNucDummy
  character(len=128) :: Filename
  character(len=80) :: VecTyp
  real(kind=wp), allocatable :: CMO(:), EDum(:), Ene(:), Occ(:)

  ! This routine is used at normal end of a RASSCF optimization, or
  ! when using the OrbOnly keyword to create orbital files.
  !---------------------------------------------------------------------
  ! NOTE:
  ! PAM 2008: Before this subroutine replaced RASREAD, the orbital energies
  ! were sent as subroutine argument when rasread was called from rasscf.
  ! But the real argument, in rasscf, was FDIAG, which turns out to be a
  ! fixed array FDIAG(mxorb) in rasscf_global.F90.
  ! Rather than checking why that array is there, and how the values are
  ! put there, I simply use the array FDIAG in rasscf_global.F90, when needing
  ! orbital energies to WrVec calls. This may need checking later...
  !--------------------------------------------------------------------*
  !     Read the JobIph file to get the required information           *
  !--------------------------------------------------------------------*
  ! PAM Jan 2014 -- do not take POTNUC from JOBIPH; take it directly
  ! from runfile, where it was stored by seward.
  iDisk = 0
  call iDaFile(JobIph,2,iToc,15,iDisk)
  iDisk = iToc(1)
  call WR_RASSCF_Info(JobIph,2,iDisk,nActEl,iSpin,nSym,stSym,nFro,nIsh,nAsh,nDel,nBas,mxSym,BName,(LenIn+8)*mxOrb,nConf,Header, &
                      144,Title,4*18*mxTit,PotNucDummy,lRoots,nRoots,iRoot,mxRoot,nRs1,nRs2,nRs3,nHole1,nElec3,iPt2,Weight)
  !if (DWSCF%do_DW) then
  !  !! recompute the weight and save on JobIph
  !  !! this is not well tested
  !  call mma_allocate(Ene,mxRoot*mxIter,Label='Ene')
  !  call Get_dArray('Last energies',Ene,lRoots)
  !  call DWSol_wgt(1,Ene,weight)
  !  call mma_deallocate(Ene)
  !  iDisk = iToc(1)
  !  call WR_RASSCF_Info(JobIph,1,iDisk,nActEl,iSpin,nSym,stSym,nFro,nIsh,nAsh,nDel,nBas,mxSym,BName,(LenIn+8)*mxOrb,nConf,Header, &
  !                      144,Title,4*18*mxTit,PotNucDummy,lRoots,nRoots,iRoot,mxRoot,nRs1,nRs2,nRs3,nHole1,nElec3,iPt2,Weight)
  !end if
  ntot = sum(nBas(:nSym))
  ntot2 = sum(nBas(:nSym)**2)
  !--------------------------------------------------------------------*
  !     Allocate CMO array                                             *
  !--------------------------------------------------------------------*
  call mma_allocate(CMO,ntot2,Label='CMO')
  call mma_allocate(Occ,ntot,Label='Occ')
  !--------------------------------------------------------------------*
  !     Make typeindex information                                     *
  !--------------------------------------------------------------------*
  if (.not. iDoGas) then
    IndType(:,:) = get_typeidx(nFro,nIsh,nRs1,nRs2,nRs3,nBas,nDel)
  else
    IndType(:,:) = get_typeidx(nFro,nIsh,nGSSH,nBas,nDel)
  end if
  !--------------------------------------------------------------------*
  !     First, write orbitals to RasOrb:
  ! IORBTYP=1 for 'Average' orbitals... Default!
  ! IORBTYP=2 for 'Canonical' orbitals.
  ! IORBTYP=3 for 'Natural' orbitals... in this case the number of roots need to be specified
  ! IORBTYP=4 for 'Spin' orbitals... in this case the number of roots need to be specified
  !--------------------------------------------------------------------*
  filename = 'RASORB'
  if (iOrbTyp /= 2) then
    iDisk = iToc(2)
    call dDaFile(JobIph,2,CMO,ntot2,iDisk)
    if (IPRLEV >= USUAL) write(u6,'(6X,3A)') 'Average orbitals are written to the ',trim(filename),' file'
    VecTyp = '* RASSCF average (pseudo-natural) orbitals'
    call dDaFile(JobIph,2,Occ,ntot,iDisk)
  else
    iDisk = iToc(9)
    call dDaFile(JobIph,2,CMO,ntot2,iDisk)
    if (IPRLEV >= USUAL) write(u6,'(6X,3A)') 'Canonical orbitals are written to the ',trim(filename),' file'
    VecTyp = '* RASSCF canonical orbitals for CASPT2'
    Occ(:) = One
  end if
  !--------------------------------------------------------------------*
  !     Write  orbitals                                                *
  !--------------------------------------------------------------------*
  LuvvVec = 50
  LuvvVec = isfreeunit(LuvvVec)
  !call WrVec(filename,LuvvVec,'COE',nSym,nBas,nBas,CMO,Occ,FDIAG,iDummy,VecTyp)
  call WrVec_(filename,LuvvVec,'COET',0,nSym,nBas,nBas,CMO,CMO,Occ,Occ,FDIAG,[E2act],indType,VecTyp,0)
  !call WrVec(filename,LuvvVec,'AI',NSYM,NBAS,NBAS,CMO,Occ,FDIAG,IndType,VecTyp)
  call WrVec_(filename,LuvvVec,'AIT',0,nSym,nBas,nBas,CMO,CMO,Occ,Occ,FDIAG,[E2act],indType,VecTyp,0)
  !--------------------------------------------------------------------*
  !     Second, write natural orbitals                                 *
  !--------------------------------------------------------------------*
  call mma_allocate(Ene,mxRoot*mxIter,Label='Ene')
  call Get_dArray('Last energies',Ene,lRoots)

  iDisk = iToc(12)
  do IRT=1,min(MAXORBOUT,LROOTS,999)
    energy = Ene(IRT)
    if (irt < 999) then
      filename = 'RASORB.'//str(IRT)
    else
      filename = 'RASORB.x'
    end if
    call dDaFile(JobIph,2,CMO,ntot2,iDisk)
    call dDaFile(JobIph,2,Occ,ntot,iDisk)
    if (IPRLEV >= USUAL) write(u6,'(6X,A,I3,3A)') 'Natural orbitals for root ',IRT,' are written to the ',trim(filename),' file'
    write(VecTyp,'(A41,I3,A3,f22.12)') '* RASSCF natural orbitals for root number',IRT,' E=',Energy
    !------------------------------------------------------------------*
    !     Write  orbitals                                              *
    !------------------------------------------------------------------*
    LuvvVec = 50
    LuvvVec = isfreeunit(LuvvVec)
    call WrVec(filename,LuvvVec,'COE',nSym,nBas,nBas,CMO,Occ,FDIAG,IndType,VecTyp)
    call WrVec(filename,LuvvVec,'AI',NSYM,NBAS,NBAS,CMO,Occ,FDIAG,IndType,VecTyp)
  end do
  call mma_deallocate(Ene)
  !--------------------------------------------------------------------*
  !     Third, write spin density orbitals                             *
  !--------------------------------------------------------------------*
  iDisk = iToc(14)
  call mma_allocate(EDum,NTot,Label='EDum')
  EDum(:) = Zero
  do IRT=1,min(MAXORBOUT,LROOTS,999)
    if (irt < 999) then
      filename = 'SPDORB.'//str(IRT)
    else
      filename = 'SPDORB.x'
    end if
    call dDaFile(JobIph,2,CMO,ntot2,iDisk)
    call dDaFile(JobIph,2,Occ,ntot,iDisk)
    if (IPRLEV >= USUAL) write(u6,'(6X,A,I3,3A)') 'Spin density orbitals for root ',irt,' are written to the ',trim(filename), &
                                                  ' file'
    write(VecTyp,'(A,I3)') '* RASSCF spin density orbitals for root number',IRT
    !------------------------------------------------------------------*
    !     Write  orbitals                                              *
    !------------------------------------------------------------------*
    LuvvVec = 50
    LuvvVec = isfreeunit(LuvvVec)
    call WrVec(filename,LuvvVec,'CEO',nSym,nBas,nBas,CMO,Occ,EDum,IndType,VecTyp)
    call WrVec(filename,LuvvVec,'AI',NSYM,NBAS,NBAS,CMO,Occ,EDum,IndType,VecTyp)
  end do
  call mma_deallocate(EDum)
  !--------------------------------------------------------------------*
  !     Normal Exit                                                    *
  !--------------------------------------------------------------------*
  call mma_deallocate(CMO)
  call mma_deallocate(Occ)

  return

end subroutine

function RAS_get_typeidx(nFro,nIsh,nRs1,nRs2,nRs3,nBas,nDel) result(typeidx)

  integer(kind=iwp) :: typeidx(7,8)
  integer(kind=iwp), intent(in) :: nFro(:), nIsh(:), nRs1(:), nRs2(:), nRs3(:), nBas(:), nDel(:)

  typeidx(1,:nSym) = nFro(:nSym)
  typeidx(2,:nSym) = nIsh(:nSym)
  typeidx(3,:nSym) = nRS1(:nSym)
  typeidx(4,:nSym) = nRS2(:nSym)
  typeidx(5,:nSym) = nRS3(:nSym)
  typeidx(7,:nSym) = nDel(:nSym)

  typeidx(6,:nSym) = 0
  typeidx(6,:nSym) = nBas(:nSym)-sum(typeidx(:,:nSym),dim=1)

end function RAS_get_typeidx

function GAS_get_typeidx(nFro,nIsh,nGSSH,nBas,nDel) result(typeidx)

  integer(kind=iwp) :: typeidx(7,8)
  integer(kind=iwp), intent(in) :: nFro(:), nIsh(:), nGSSH(:,:), nBas(:), nDel(:)

  typeidx(1,:nSym) = nFro(:nSym)
  typeidx(2,:nSym) = nIsh(:nSym)
  typeidx(3,:nSym) = 0
  typeidx(4,:nSym) = sum(nGssh(1:nGAS,:nSym),dim=1)
  typeidx(5,:nSym) = 0
  typeidx(7,:nSym) = nDel(:nSym)

  typeidx(6,:nSym) = 0
  typeidx(6,:nSym) = nBas(:nSym)-sum(typeidx(:,:nSym),dim=1)

end function GAS_get_typeidx

subroutine putOrbFile(CMO,orbital_E,iDoGAS)

  real(kind=wp), intent(in) :: CMO(:), orbital_E(:)
  logical(kind=iwp), intent(in) :: iDoGAS
  integer(kind=iwp) :: file_id, typeidx(7,8)
  real(kind=wp), allocatable :: occ_number(:)
  character(len=*), parameter :: filename = 'ORTHORB', orbfile_title = 'Orbitals after Orthonormalization.'
  integer(kind=iwp), parameter :: arbitrary_magic_number = 50

  file_id = isfreeunit(arbitrary_magic_number)
  if (.not. iDoGas) then
    typeidx = get_typeidx(nFro,nIsh,nRs1,nRs2,nRs3,nBas,nDel)

  else
    typeidx = get_typeidx(nFro,nIsh,nGSSH,nBas,nDel)
  end if

  ! TODO(Oskar): Implement proper occupation number reading.
  call mma_allocate(occ_number,nTot)
  occ_number(:) = One
  call WrVec(filename,file_id,'COIE',nSym,nBas,nBas,CMO,occ_number,orbital_E,typeidx,orbfile_title)
  call mma_deallocate(occ_number)

end subroutine putOrbFile

end module write_orbital_files
