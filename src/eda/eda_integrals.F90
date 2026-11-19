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

module EDA_Integrals

use Definitions, only: wp, iwp
use Index_Functions, only: nTri_Elem
use Constants, only: Zero
use Para_Info, only: nProcs, MyRank, Is_Real_Par

implicit none
private
public :: EDA_Read_OneElectron, EDA_TwoElectron_Energies

contains

subroutine EDA_Read_OneElectron(nBas,H,S,Enuc,ierr,Message)

  use OneDat, only: sNoNuc, sNoOri

  integer(kind=iwp), intent(in) :: nBas
  real(kind=wp), intent(out) :: H(nBas,nBas), S(nBas,nBas)
  real(kind=wp), intent(out) :: Enuc
  integer(kind=iwp), intent(out) :: ierr
  character(len=*), intent(out) :: Message

  integer(kind=iwp) :: nPair, rc, CloseRC, iOpt, iComp, iSym, LuOne, AnyError
  logical :: Opened, Failed, IsParallel
  character(len=8) :: Label
  real(kind=wp), allocatable :: Buf(:)

  IsParallel = (nProcs > 1) .and. Is_Real_Par()

  ierr = 0
  Message = ''
  H(:,:) = Zero
  S(:,:) = Zero
  Enuc = Zero
  Opened = .false.
  rc = 0
  CloseRC = 0
  LuOne = 41
  nPair = nTri_Elem(nBas)

  ! ONEINT access is local in the Molcas reader layer.  Read the complete
  ! one-electron data only on rank 0 and broadcast it by a global sum from
  ! zero-filled arrays on the remaining ranks.
  if ((.not. IsParallel .or. MyRank == 0)) then
    allocate(Buf(nPair+4))

    call OpnOne(rc,0,'ONEINT',LuOne)
    if (rc == 0) then
      Opened = .true.
      iOpt = ibset(ibset(0,sNoOri),sNoNuc)
      iComp = 1
      iSym = 1

      Label = 'OneHam  '
      call RdOne(rc,iOpt,Label,iComp,Buf,iSym)
      if (rc == 0) call Square(Buf,H,nBas,1,nBas)

      if (rc == 0) then
        Label = 'Mltpl  0'
        call RdOne(rc,iOpt,Label,iComp,Buf,iSym)
        if (rc == 0) call Square(Buf,S,nBas,1,nBas)
      end if

      if (rc == 0) call Get_dScalar('PotNuc',Enuc)
    end if

    if (Opened) then
      call ClsOne(CloseRC,0)
      Opened = .false.
      if ((rc == 0) .and. (CloseRC /= 0)) rc = CloseRC
    end if
    if (allocated(Buf)) deallocate(Buf)
  end if

  Failed = (rc /= 0)
  AnyError = 0
  if (Failed) AnyError = 1
  if (IsParallel) call gaIgOP_SCAL(AnyError,'max')
  if (AnyError /= 0) then
    ierr = 1
    if (Failed) then
      Message = 'Unable to read ONEINT for EDA.'
    else
      Message = 'ONEINT reading failed on another MPI rank.'
    end if
    return
  end if

  if (IsParallel) then
    call GADGOp(H(1,1),nBas*nBas,'+')
    call GADGOp(S(1,1),nBas*nBas,'+')
    call GADGOp_SCAL(Enuc,'+')
  end if

end subroutine EDA_Read_OneElectron

subroutine EDA_TwoElectron_Energies(nBas,nFrag,Praw,Porth,Pfrag, &
                                     Jraw,Kraw,Jorth,Korth,Kfrag,ierr,Message)

  integer(kind=iwp), intent(in) :: nBas, nFrag
  real(kind=wp), intent(in) :: Praw(nBas,nBas), Porth(nBas,nBas)
  real(kind=wp), intent(in) :: Pfrag(nBas,nBas,nFrag)
  real(kind=wp), intent(out) :: Jraw,Kraw,Jorth,Korth
  real(kind=wp), intent(out) :: Kfrag(nFrag)
  integer(kind=iwp), intent(out) :: ierr
  character(len=*), intent(out) :: Message

  integer(kind=iwp) :: nPair, rc, CloseRC, LuOrd, AnyError, iPass
  integer(kind=iwp) :: FileNSym, FileBas(8), FileSkip(8)
  integer(kind=iwp) :: i,j,k,l,q,nMat,iOpt,iFrag
  integer(kind=iwp) :: mij,mkl
  logical(kind=iwp) :: Exists, Square
  logical :: Opened, Failed, IsParallel
  real(kind=wp) :: v, JFactor
  real(kind=wp), allocatable :: Buf(:), Totals(:)

  IsParallel = (nProcs > 1) .and. Is_Real_Par()

  ierr = 1
  Message = 'Unable to read conventional two-electron integrals for EDA.'
  Jraw = Zero
  Kraw = Zero
  Jorth = Zero
  Korth = Zero
  Kfrag(:) = Zero
  Opened = .false.
  rc = 0
  CloseRC = 0
  LuOrd = 43
  nPair = nTri_Elem(nBas)

  call f_Inquire('ORDINT',Exists)
  AnyError = 0
  if (.not. logical(Exists)) AnyError = 1
  if (IsParallel) call gaIgOP_SCAL(AnyError,'max')
  if (AnyError /= 0) then
    Message = 'EDA requires stored ORDINT integrals: use &SEWARD with NoCholesky.'
    return
  end if

  allocate(Buf(nPair+1))
  reading: do iPass=1,1

  call OpnOrd(rc,0,'ORDINT',LuOrd)
  if (rc == 0) Opened = .true.
  AnyError = 0
  if (rc /= 0) AnyError = 1
  if (IsParallel) call gaIgOP_SCAL(AnyError,'max')
  if (AnyError /= 0) then
    Message = 'Unable to open conventional ORDINT integrals for EDA.'
    exit reading
  end if

  call GetOrd(rc,Square,FileNSym,FileBas,FileSkip)
  Failed = (rc /= 0)
  if (.not. Failed) Failed = (FileNSym /= 1) .or. (FileBas(1) /= nBas) .or. (FileSkip(1) /= 0)
  AnyError = 0
  if (Failed) AnyError = 1
  if (IsParallel) call gaIgOP_SCAL(AnyError,'max')
  if (AnyError /= 0) then
    Message = 'ORDINT dimensions do not match the C1 EDA supermolecule reference.'
    exit reading
  end if

  ! In a real MPI run ORDINT is distributed over the ranks.  Unlike an
  ! AO->MO transformation, each EDA J/K quantity is a linear functional
  ! of the AO integrals.  It is therefore sufficient for every rank to
  ! contract its local ORDINT contribution and globally sum only the final
  ! scalar energies.  No reconstruction of complete AO rows is required.
  iOpt = 1
  do i=1,nBas
    do j=1,i
      Buf(:) = Zero
      rc = 0
      nMat = 0
      call RdOrd_(rc,iOpt,1,1,1,1,Buf,nPair+1,nMat)
      iOpt = 2
      Failed = (rc /= 0) .or. (nMat /= 1)
      AnyError = 0
      if (Failed) AnyError = 1
      if (IsParallel) call gaIgOP_SCAL(AnyError,'max')
      if (AnyError /= 0) then
        Message = 'Failed to read a complete conventional ORDINT pair row for EDA.'
        exit reading
      end if

      if (i == j) then
        mij = 1
      else
        mij = 2
      end if

      q = 0
      do k=1,nBas
        do l=1,k
          q = q+1
          v = Buf(q)
          if (v == Zero) cycle

          if (k == l) then
            mkl = 1
          else
            mkl = 2
          end if
          JFactor = real(mij*mkl,wp)*v

          Jraw = Jraw + JFactor*Praw(i,j)*Praw(k,l)
          Jorth = Jorth + JFactor*Porth(i,j)*Porth(k,l)

          Kraw = Kraw + v*ExchangeMultiplicity(Praw,i,j,k,l)
          Korth = Korth + v*ExchangeMultiplicity(Porth,i,j,k,l)
          do iFrag=1,nFrag
            Kfrag(iFrag) = Kfrag(iFrag) + v*ExchangeMultiplicity(Pfrag(:,:,iFrag),i,j,k,l)
          end do
        end do
      end do
    end do
  end do

  rc = 0
  if (Opened) then
    call ClsOrd(CloseRC)
    Opened = .false.
    if (CloseRC /= 0) rc = CloseRC
  end if
  AnyError = 0
  if (rc /= 0) AnyError = 1
  if (IsParallel) call gaIgOP_SCAL(AnyError,'max')
  if (AnyError /= 0) then
    Message = 'Unable to close ORDINT after EDA.'
    exit reading
  end if

  if (IsParallel) then
    allocate(Totals(4+nFrag))
    Totals(1) = Jraw
    Totals(2) = Kraw
    Totals(3) = Jorth
    Totals(4) = Korth
    Totals(5:4+nFrag) = Kfrag(1:nFrag)
    call GADGOp(Totals(1),4+nFrag,'+')
    Jraw = Totals(1)
    Kraw = Totals(2)
    Jorth = Totals(3)
    Korth = Totals(4)
    Kfrag(1:nFrag) = Totals(5:4+nFrag)
    deallocate(Totals)
  end if

  ierr = 0
  Message = ''
  end do reading

  CloseRC = 0
  if (Opened) then
    call ClsOrd(CloseRC)
    Opened = .false.
  end if
  AnyError = 0
  if (CloseRC /= 0) AnyError = 1
  if (IsParallel) call gaIgOP_SCAL(AnyError,'max')
  if (AnyError /= 0) then
    ierr = 1
    if (len_trim(Message) == 0) Message = 'Unable to close ORDINT after EDA.'
  end if
  if (allocated(Buf)) deallocate(Buf)
  if (allocated(Totals)) deallocate(Totals)

contains

  real(kind=wp) function ExchangeMultiplicity(P,a,b,c,d)

    real(kind=wp), intent(in) :: P(nBas,nBas)
    integer(kind=iwp), intent(in) :: a,b,c,d

    ExchangeMultiplicity = P(a,c)*P(b,d)
    if (a /= b) ExchangeMultiplicity = ExchangeMultiplicity + P(b,c)*P(a,d)
    if (c /= d) ExchangeMultiplicity = ExchangeMultiplicity + P(a,d)*P(b,c)
    if ((a /= b) .and. (c /= d)) ExchangeMultiplicity = ExchangeMultiplicity + P(b,d)*P(a,c)

  end function ExchangeMultiplicity

end subroutine EDA_TwoElectron_Energies

end module EDA_Integrals
