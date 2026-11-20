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
! Copyright (C) 2026, Okan Koeksal                                     *
!***********************************************************************

module EDA_Energy

use Definitions, only: wp, iwp
use Constants, only: Zero, One, Two
use EDA_Reference, only: EDA_Reference_Type
use Para_Info, only: nProcs, Is_Real_Par
use EDA_Integrals, only: EDA_Read_OneElectron, EDA_TwoElectron_Energies

implicit none
private

public :: EDA_Results_Type, EDA_Compute_ClosedShell

! P. Su and H. Li, J. Chem. Phys. 131, 014102 (2009),
! DOI: 10.1063/1.3159673.
! Paper notation used below:
!   EX1 = E_X^(1), EX2 = E_X^(2), EX3 = E_X^(3)
!   DeltaEele, DeltaEex, DeltaErep, DeltaEpol, DeltaEHF, DeltaEdisp.
type :: EDA_Results_Type
  real(kind=wp) :: EX1 = Zero
  real(kind=wp) :: EX2 = Zero
  real(kind=wp) :: EX3 = Zero
  real(kind=wp) :: DeltaEele = Zero
  real(kind=wp) :: DeltaEex = Zero
  real(kind=wp) :: DeltaErep = Zero
  real(kind=wp) :: DeltaEpol = Zero
  real(kind=wp) :: DeltaEHF = Zero
  real(kind=wp) :: DeltaEdisp = Zero
  real(kind=wp) :: DeltaEMP2 = Zero
  real(kind=wp) :: Closure = Zero
  logical :: HasMP2 = .false.
end type EDA_Results_Type

contains

subroutine EDA_Compute_ClosedShell(Super,Frag,nFrag,FragBasis,Results,ierr,Message)

  type(EDA_Reference_Type), intent(in) :: Super
  type(EDA_Reference_Type), intent(in) :: Frag(:)
  integer(kind=iwp), intent(in) :: nFrag
  character(len=*), intent(in) :: FragBasis
  type(EDA_Results_Type), intent(out) :: Results
  integer(kind=iwp), intent(out) :: ierr
  character(len=*), intent(out) :: Message

  integer(kind=iwp) :: nBas,nOcc,nOccSum,nBasSum
  integer(kind=iwp) :: iFrag, i0, i1, mu, nu, p
  integer(kind=iwp) :: info, lWork, iPass
  real(kind=wp) :: Enuc,OneRaw,OneOrth,Jraw,Kraw,Jorth,Korth
  real(kind=wp) :: SumEHF,SumEMP2,OverlapMin
  real(kind=wp), allocatable :: H(:,:),S(:,:),C(:,:),Tmp(:,:),OccS(:,:),Eig(:),Work(:),InvOccS(:,:)
  real(kind=wp), allocatable :: Praw(:,:),Porth(:,:),Pfrag(:,:,:),Kfrag(:)
  logical :: AllMP2, IsParallel
  integer(kind=iwp) :: AnyError

  ierr = 1
  Message = ''
  Results = EDA_Results_Type()
  IsParallel = (nProcs > 1) .and. Is_Real_Par()

  nBas = Super%nBas
  nOcc = Super%nOcc
  AnyError = 0
  if ((nBas < 1) .or. (nOcc < 1)) AnyError = 1
  if (IsParallel) call gaIgOP_SCAL(AnyError,'max')
  if (AnyError /= 0) then
    Message = 'Invalid supermolecule RHF dimensions.'
    return
  end if

  nOccSum = 0
  nBasSum = 0
  do iFrag=1,nFrag
    nOccSum = nOccSum + Frag(iFrag)%nOcc
    nBasSum = nBasSum + Frag(iFrag)%nBas
  end do
  AnyError = 0
  if (nOccSum /= nOcc) AnyError = 1
  if (IsParallel) call gaIgOP_SCAL(AnyError,'max')
  if (AnyError /= 0) then
    Message = 'Sum of fragment occupied orbitals differs from the supermolecule RHF occupation.'
    return
  end if

  if (trim(FragBasis) == 'OWN') then
    AnyError = 0
    if (nBasSum /= nBas) AnyError = 1
    if (IsParallel) call gaIgOP_SCAL(AnyError,'max')
    if (AnyError /= 0) then
      Message = 'FBASIS=OWN requires the sum of fragment AO dimensions to equal the supermolecule AO dimension.'
      return
    end if
  else if (trim(FragBasis) == 'SUPER') then
    do iFrag=1,nFrag
      AnyError = 0
      if (Frag(iFrag)%nBas /= nBas) AnyError = 1
      if (IsParallel) call gaIgOP_SCAL(AnyError,'max')
      if (AnyError /= 0) then
        Message = 'FBASIS=SUPER requires every fragment RunFile to use the full supermolecule AO basis.'
        return
      end if
    end do
  else
    Message = 'Internal EDA error: unknown fragment-basis convention.'
    return
  end if

  allocate(H(nBas,nBas),S(nBas,nBas),C(nBas,nOcc),Tmp(nBas,nOcc))
  allocate(OccS(nOcc,nOcc),Eig(nOcc),InvOccS(nOcc,nOcc))
  allocate(Praw(nBas,nBas),Porth(nBas,nBas),Pfrag(nBas,nBas,nFrag),Kfrag(nFrag))

  processing: do iPass=1,1

  call EDA_Read_OneElectron(nBas,H,S,Enuc,ierr,Message)
  if (ierr /= 0) exit processing

  C(:,:) = Zero
  Pfrag(:,:,:) = Zero
  i0 = 1
  i1 = 0

  if (trim(FragBasis) == 'OWN') then
    p = 0
    do iFrag=1,nFrag
      i1 = i0 + Frag(iFrag)%nOcc - 1
      do nu=1,Frag(iFrag)%nOcc
        do mu=1,Frag(iFrag)%nBas
          C(p+mu,i0+nu-1) = Frag(iFrag)%COcc(mu,nu)
        end do
      end do
      call DGEMM_('N','T',nBas,nBas,Frag(iFrag)%nOcc,Two, &
                  C(1,i0),nBas,C(1,i0),nBas,Zero,Pfrag(1,1,iFrag),nBas)
      p = p + Frag(iFrag)%nBas
      i0 = i1 + 1
    end do
  else
    do iFrag=1,nFrag
      i1 = i0 + Frag(iFrag)%nOcc - 1
      C(:,i0:i1) = Frag(iFrag)%COcc(:,:)
      call DGEMM_('N','T',nBas,nBas,Frag(iFrag)%nOcc,Two, &
                  C(1,i0),nBas,C(1,i0),nBas,Zero,Pfrag(1,1,iFrag),nBas)
      i0 = i1 + 1
    end do
  end if

  Praw(:,:) = Zero
  do iFrag=1,nFrag
    Praw(:,:) = Praw(:,:) + Pfrag(:,:,iFrag)
  end do

  ! Occupied-orbital overlap: O = C^T S C.
  call DGEMM_('N','N',nBas,nOcc,nBas,One,S,nBas,C,nBas,Zero,Tmp,nBas)
  call DGEMM_('T','N',nOcc,nOcc,nBas,One,C,nBas,Tmp,nBas,Zero,OccS,nOcc)

  ! Symmetric eigendecomposition O = U diag(s) U^T, then O^-1.
  allocate(Work(1))
  lWork = -1
  call DSYEV_('V','U',nOcc,OccS,nOcc,Eig,Work,lWork,info)
  AnyError = 0
  if (info /= 0) AnyError = 1
  if (IsParallel) call gaIgOP_SCAL(AnyError,'max')
  if (AnyError /= 0) then
    Message = 'Failed to query the occupied-overlap diagonalization workspace.'
    exit processing
  end if
  lWork = max(1_iwp,int(Work(1),kind=iwp))
  deallocate(Work)
  allocate(Work(lWork))
  call DSYEV_('V','U',nOcc,OccS,nOcc,Eig,Work,lWork,info)
  AnyError = 0
  if (info /= 0) AnyError = 1
  if (IsParallel) call gaIgOP_SCAL(AnyError,'max')
  if (AnyError /= 0) then
    Message = 'Failed to diagonalize the occupied fragment-orbital overlap matrix.'
    exit processing
  end if

  OverlapMin = minval(Eig)
  AnyError = 0
  if (OverlapMin <= 1.0e-10_wp) AnyError = 1
  if (IsParallel) call gaIgOP_SCAL(AnyError,'max')
  if (AnyError /= 0) then
    Message = 'Fragment occupied orbitals are linearly dependent in the supermolecule overlap metric.'
    exit processing
  end if

  InvOccS(:,:) = Zero
  do p=1,nOcc
    do nu=1,nOcc
      do mu=1,nOcc
        InvOccS(mu,nu) = InvOccS(mu,nu) + OccS(mu,p)*OccS(nu,p)/Eig(p)
      end do
    end do
  end do

  ! Porth = 2 C O^-1 C^T.
  call DGEMM_('N','N',nBas,nOcc,nOcc,One,C,nBas,InvOccS,nOcc,Zero,Tmp,nBas)
  call DGEMM_('N','T',nBas,nBas,nOcc,Two,Tmp,nBas,C,nBas,Zero,Porth,nBas)

  OneRaw = sum(Praw*H)
  OneOrth = sum(Porth*H)

  call EDA_TwoElectron_Energies(nBas,nFrag,Praw,Porth,Pfrag,Jraw,Kraw,Jorth,Korth,Kfrag,ierr,Message)
  if (ierr /= 0) exit processing

  Results%EX1 = OneRaw + 0.5_wp*Jraw - 0.25_wp*sum(Kfrag) + Enuc
  Results%EX2 = OneRaw + 0.5_wp*Jraw - 0.25_wp*Kraw + Enuc
  Results%EX3 = OneOrth + 0.5_wp*Jorth - 0.25_wp*Korth + Enuc

  SumEHF = Zero
  SumEMP2 = Zero
  AllMP2 = Super%HasMP2
  do iFrag=1,nFrag
    SumEHF = SumEHF + Frag(iFrag)%EHF
    SumEMP2 = SumEMP2 + Frag(iFrag)%ETotal
    AllMP2 = AllMP2 .and. Frag(iFrag)%HasMP2
  end do

  Results%DeltaEele = Results%EX1 - SumEHF
  Results%DeltaEex = Results%EX2 - Results%EX1
  Results%DeltaErep = Results%EX3 - Results%EX2
  Results%DeltaEpol = Super%EHF - Results%EX3
  Results%DeltaEHF = Super%EHF - SumEHF
  Results%Closure = Results%DeltaEHF - (Results%DeltaEele + Results%DeltaEex + &
                    Results%DeltaErep + Results%DeltaEpol)

  Results%HasMP2 = AllMP2
  if (AllMP2) then
    Results%DeltaEMP2 = Super%ETotal - SumEMP2
    Results%DeltaEdisp = Results%DeltaEMP2 - Results%DeltaEHF
  end if

  ierr = 0
  Message = ''
  end do processing

  if (allocated(H)) deallocate(H)
  if (allocated(S)) deallocate(S)
  if (allocated(C)) deallocate(C)
  if (allocated(Tmp)) deallocate(Tmp)
  if (allocated(OccS)) deallocate(OccS)
  if (allocated(Eig)) deallocate(Eig)
  if (allocated(Work)) deallocate(Work)
  if (allocated(InvOccS)) deallocate(InvOccS)
  if (allocated(Praw)) deallocate(Praw)
  if (allocated(Porth)) deallocate(Porth)
  if (allocated(Pfrag)) deallocate(Pfrag)
  if (allocated(Kfrag)) deallocate(Kfrag)

end subroutine EDA_Compute_ClosedShell

end module EDA_Energy
