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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!               2003, Valera Veryazov                                  *
!               1998,2022, Roland Lindh                                *
!***********************************************************************

subroutine SOrb(LuOrb,SIntTh,iTerm)
!***********************************************************************
!                                                                      *
!     purpose: Get starting orbitals from:                             *
!              -1) default choice                                      *
!               0) diagonalizaton of the core                          *
!               1) via intermediate calculation of HF AOs              *
!               2) input orbitals                                      *
!               3) input density matrix                                *
!                                                                      *
!***********************************************************************

use InfSCF, only: CMO, DoCholesky, EOrb, FockAO, InVec, KSDFT, nBas, nBB, nBB, nBT, nD, nnB, nOrb, nSym, OccNo, One_Grid, OneHam, &
                  Ovrlp, SCF_FileOrb, ScrFac, Scrmbl, StVec, TrM
#ifdef _HDF5_
use mh5, only: mh5_close_file
use InfSCF, only: FileOrb_ID, IsHDF5
#endif
use SCFFiles, only: LuOut
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: LuOrb
real(kind=wp), intent(inout) :: SIntTh
integer(kind=iwp), intent(out) :: iTerm
integer(kind=iwp) :: iD, IsUHF, nData
logical(kind=iwp) :: found, FstItr, Retry
character(len=512) :: FName
character(len=80) :: KSDFT_save

call DecideonCholesky(DoCholesky)
! Cholesky and NDDO are incompatible
if (DoCholesky .and. (InVec == 1)) then
  call WarningMessage(1,' In SORB: Cholesky and NDDO not implemented !!!;  NDDO option ignored')
  InVec = -1
end if

! Is default clause chosen?

if (InVec == -1) then
  call qpg_darray('SCF orbitals',found,ndata)
  if (found .and. (nData == nBB)) then
    call qpg_darray('OrbE',found,ndata)
    if (found) InVec = 8
  end if
end if
if (InVec == -1) then
  call qpg_darray('Guessorb',found,ndata)
  if (found .and. (nData == nBB)) then
    call qpg_darray('Guessorb energies',found,ndata)
    if (found) InVec = 9
  end if
end if
if (InVec == -1) InVec = 0
!                                                                      *
!***********************************************************************
!                                                                      *
! Has the user selected a method?

do
  Retry = .false.
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  select case (InVec)

    case (0)

      ! Diagonalize core
      call Start0()
      !                                                                *
      !*****************************************************************
      !                                                                *
    case (1)

      ! HF AO orbitals as intermediate step...

      ! NDDO, always none-DFT

      call SwiOpt(.false.,OneHam,Ovrlp,nBT,CMO,nBB,nD)
      call Start0()
      InVec = 0
      call SOrbCHk(OneHam,FockAO,nBT,nD)
      KSDFT_save = KSDFT
      KSDFT = 'SCF'
      call WrInp_SCF(SIntTh)
      FstItr = .true.
      call WfCtl_SCF(iTerm,'NDDO      ',FstItr,SIntTh)
      KSDFT = KSDFT_save
      call Free_TLists()
      if (iTerm /= 0) call Quit(iTerm)
      write(u6,*)
      write(u6,'(A)') 'Generation of NDDO vectors completed!'
      write(u6,*)
      write(u6,*) '2nd step: optimizing HF MOs...'
      write(u6,*) '------------------------------'
      call SwiOpt(.true.,OneHam,Ovrlp,nBT,CMO,nBB,nD)
      ! Reset to to start from the current MO set
      call Init_SCF()
      InVec = 5
      !IFG: I presume the arguments after LuOut in these two calls are correct,
      !     they were missing!
      if (nD == 1) then
        FName = 'SCFORB'
        call Start2(FName,LuOut,CMO,nBB,nD,Ovrlp,nBT,EOrb,OccNo,nnB)
      else
        FName = 'UHFORB'
        call Start2(FName,LuOut,CMO,nBB,nD,Ovrlp,nBT,EOrb,OccNo,nnB)
      end if
      !                                                                *
      !*****************************************************************
      !                                                                *
    case (2)

      ! Read INPORB
      One_Grid = .true.
      FName = SCF_FileOrb
      call Start2(FName,LuOrb,CMO,nBB,nD,Ovrlp,nBT,EOrb,OccNo,nnB)
      !                                                                *
      !*****************************************************************
      !                                                                *
    case (3)

      ! Read COMOLD
      One_Grid = .true.
      call Start3(CMO,TrM,nBB,nD,OneHam,Ovrlp,nBT)
      !                                                                *
      !*****************************************************************
      !                                                                *
    case (6)

      write(u6,*)
      write(u6,*) '     Constrained SCF calculation '
      write(u6,*)
      StVec = 'Constrained orbitals'
      One_Grid = .true.
      FName = SCF_FileOrb
      call Chk_Vec_UHF(FName,LuOrb,isUHF)
      if (IsUHF == 1) then
        InVec = 2
        Retry = .true.
      end if
      call Start6(FName,LuOrb,CMO,nBB,nD,EOrb,OccNo,nnB)
      !                                                                *
      !*****************************************************************
      !                                                                *
    case (8)

      StVec = 'Detected old SCF orbitals'
      One_Grid = .true.
      call start0y(CMO,nBB,nD,EOrb,nnB)
      !                                                                *
      !*****************************************************************
      !                                                                *
    case (9)

      StVec = 'Detected guessorb starting orbitals'
      !One_Grid = .true.
      call start0x(CMO,nBB,nD,EOrb,nnB)
      !                                                                *
      !*****************************************************************
      !                                                                *
    case default

      write(u6,*) 'Illegal inVec value:',InVec
      call Abend()
      !                                                                *
      !*****************************************************************
      !                                                                *
  end select
  if (.not. Retry) exit
end do
!                                                                      *
!***********************************************************************
!                                                                      *
if (Scrmbl) then
  do iD=1,nD
    call Scram(CMO(1,iD),nSym,nBas,nOrb,ScrFac)
  end do
end if

call SOrbCHk(OneHam,FockAO,nBT,nD)
#ifdef _HDF5_
if (isHDF5) call mh5_close_file(fileorb_id)
#endif

contains

subroutine Start0()
!***********************************************************************
!                                                                      *
!     purpose: Get starting orbitals from diagonalization of the core. *
!              The result is stored in CMO. Those orbitals are         *
!              then optimized in SubRoutine WfCtl. A set of orthogonal,*
!              symmetry adapted AOs is stored in TrM.                  *
!                                                                      *
!***********************************************************************

  use InfSCF, only: nBO, nOcc, nnB

  integer(kind=iwp) :: iD

  !--------------------------------------------------------------------*
  !     Start                                                          *
  !--------------------------------------------------------------------*

  ! Form transformation matrix
  call TrGen(TrM(1,1),nBB,Ovrlp,OneHam,nBT)
  if (nD == 2) TrM(:,2) = TrM(:,1)

  ! Diagonalize core
  do iD=1,nD
    call DCore(OneHam,nBT,CMO(1,iD),TrM(1,iD),nBO,EOrb(1,iD),nnB,nOcc(1,iD),Ovrlp)
  end do

  !--------------------------------------------------------------------*
  !     Exit                                                           *
  !--------------------------------------------------------------------*
  return

end subroutine Start0

subroutine Start0x(CMO,mBB,nD,E_Or,mmB)
!***********************************************************************
!                                                                      *
! This routine reads start orbitals generated by guessorb.             *
!                                                                      *
!***********************************************************************

  use InfSCF, only: nDel

  integer(kind=iwp) :: mBB, nD, mmB
  real(kind=wp) :: CMO(mBB,nD), E_Or(mmB,nD)
  integer(kind=iwp) :: iD, iRC, nSum

  !--------------------------------------------------------------------*
  ! Get start orbitals.                                                *
  !--------------------------------------------------------------------*

  iRC = 0
  call qpg_darray('Guessorb',Found,ndata)
  if (Found) then
    if (nData /= mBB) then
      write(u6,*) 'Start0x: nData /= mBB'
      write(u6,*) '         nData=',nData
      write(u6,*) '         mBB  =',mBB
      call Abend()
    end if
    call get_darray('Guessorb',CMO(:,1),ndata)
  else
    iRC = -1
  end if
  if (iRC /= 0) then
    write(u6,*) 'Start0x: no orbitals found!'
    call Abend()
  end if

  iRC = 0
  call qpg_darray('Guessorb energies',found,ndata)
  if (Found) then
    if (nData /= mmB) then
      write(u6,*) 'Start0x: nData /= mmB'
      write(u6,*) '         nData=',nData
      write(u6,*) '         mmB  =',mmB
      call Abend()
    end if
    call get_darray('Guessorb energies',E_Or(:,1),ndata)
  else
    iRC = -1
  end if
  if (iRC /= 0) then
    write(u6,*) 'Start0x: no energies found!'
    call Abend()
  end if

  if (nD == 2) then
    CMO(:,2) = CMO(:,1)
    E_Or(:,2) = E_Or(:,1)
  end if

  call qpg_iarray('nDel_go',Found,ndata)
  nSum = 0
  if (Found) then
    call Get_iArray('nDel_go',nDel,ndata)
    call Put_iArray('nDel',nDel,ndata)
    nSum = sum(nDel(1:nSym))
  end if

  if (nSum > 0) then
    nOrb(1:nSym) = nBas(1:nSym)-nDel(1:nSym)
    do iD=1,nD
      call TrimCMO(CMO(:,iD),nSym,nBas,nOrb)
      call TrimEor(E_Or(:,iD),nSym,nBas,nOrb)
    end do
  end if

  !--------------------------------------------------------------------*
  ! Done                                                               *
  !--------------------------------------------------------------------*
  return

end subroutine Start0x

subroutine Start0y(CMO,mBB,nD,E_Or,mmB)
!***********************************************************************
!                                                                      *
! This routine reads old SCf orbitals as start orbitals.               *
!                                                                      *
!***********************************************************************

  use InfSCF, only: nDel

  integer(kind=iwp) :: mBB, nD, mmB
  real(kind=wp) :: CMO(mBB,nD), E_Or(mmB,nD)
  integer(kind=iwp) :: iD, nSum

  !--------------------------------------------------------------------*
  ! Get start orbitals.                                                *
  !--------------------------------------------------------------------*

  call qpg_darray('SCF orbitals',found,ndata)
  if (Found) call get_darray('SCF orbitals',CMO(:,1),ndata)
  call qpg_darray('OrbE',found,ndata)
  if (Found) call get_darray('OrbE',E_Or(:,1),ndata)

  if (nD == 2) then

    CMO(:,2) = CMO(:,1)
    E_Or(:,2) = E_Or(:,1)

    call qpg_darray('SCF orbitals_ab',found,ndata)
    if (Found) call get_darray('SCF orbitals_ab',CMO(:,2),ndata)
    call qpg_darray('OrbE_ab',found,ndata)
    if (found) call get_darray('OrbE_ab',E_Or(:,2),ndata)
  end if

  call qpg_iarray('nDel',Found,ndata)
  nSum = 0
  if (Found) then
    call Get_iArray('nDel',nDel,ndata)
    nSum = sum(nDel(1:nSym))
  end if
  if (nSum > 0) then
    nOrb(1:nSym) = nBas(1:nSym)-nDel(1:nSym)
    do iD=1,nD
      call TrimCMO(CMO(:,iD),nSym,nBas,nOrb)
      call TrimEor(E_or(:,iD),nSym,nBas,nOrb)
    end do
  end if
  !--------------------------------------------------------------------*
  ! Done                                                               *
  !--------------------------------------------------------------------*
  return

end subroutine Start0y

subroutine SOrbChk(OneHam,Fock,mBT,nD)

  integer(kind=iwp) :: mBT, nD
  real(kind=wp) :: OneHam(mBT), Fock(mBT,nD)
  integer(kind=iwp) :: iD
  real(kind=wp) :: Whatever

  do iD=1,nD
    ! Check orthonormality of start orbitals
    call ChkOrt(iD,Whatever)

    ! Form the first Fock matrix
    Fock(:,iD) = OneHam(:)
  end do

  return

end subroutine SOrbChk

end subroutine SOrb
