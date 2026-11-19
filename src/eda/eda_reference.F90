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

module EDA_Reference

use Definitions, only: wp, iwp
use Constants, only: Zero

implicit none
private

type, public :: EDA_Reference_Type
  integer(kind=iwp) :: nBas = 0
  integer(kind=iwp) :: nOrb = 0
  integer(kind=iwp) :: nOcc = 0
  integer(kind=iwp) :: nAtoms = 0
  real(kind=wp) :: EHF = Zero
  real(kind=wp) :: ETotal = Zero
  logical :: HasMP2 = .false.
  character(len=8) :: Method = '        '
  real(kind=wp), allocatable :: COcc(:,:)
end type EDA_Reference_Type

public :: EDA_Read_Super, EDA_Read_Fragment, EDA_Free_Reference

contains

subroutine EDA_Read_Super(Ref,ierr,Message)

  type(EDA_Reference_Type), intent(inout) :: Ref
  integer(kind=iwp), intent(out) :: ierr
  character(len=*), intent(out) :: Message

  call Read_Current_RunFile(Ref,ierr,Message)

end subroutine EDA_Read_Super

subroutine EDA_Read_Fragment(iFrag,Ref,ierr,Message)

  integer(kind=iwp), intent(in) :: iFrag
  type(EDA_Reference_Type), intent(inout) :: Ref
  integer(kind=iwp), intent(out) :: ierr
  character(len=*), intent(out) :: Message
  logical(kind=iwp) :: Exists
  character(len=8) :: RunName

  ierr = 1
  Message = ''

  if (iFrag < 1) then
    Message = 'EDA fragment index must be positive.'
    return
  end if

  call Fragment_RunName(iFrag,RunName,ierr)
  if (ierr /= 0) then
    Message = 'Unable to encode the EDA fragment RunFile name.'
    return
  end if

  call f_Inquire(trim(RunName),Exists)
  if (.not. Exists) then
    write(Message,'(A,I0,2A)') 'Missing EDA fragment RunFile ',iFrag,' (alias ',trim(RunName)//').'
    return
  end if

  call NameRun(trim(RunName))
  call Read_Current_RunFile(Ref,ierr,Message)
  call NameRun('#Pop')

end subroutine EDA_Read_Fragment

subroutine Fragment_RunName(iFrag,RunName,ierr)

  integer(kind=iwp), intent(in) :: iFrag
  character(len=8), intent(out) :: RunName
  integer(kind=iwp), intent(out) :: ierr
  integer(kind=iwp) :: Value
  integer :: Pos, Digit
  character(len=*), parameter :: Digits = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  RunName = '        '
  ierr = 0

  ! Preserve the historical aliases for ordinary fragment counts.
  if (iFrag <= 9) then
    write(RunName,'("EDAFRG",I1)') iFrag
    return
  else if (iFrag <= 999) then
    write(RunName,'("EDAFR",I0)') iFrag
    return
  end if

  ! For larger counts, encode the positive fragment index in base 36.
  ! Six base-36 digits with the two-character prefix "EF" cover more
  ! than two billion fragment indices while keeping an eight-character
  ! RunFile alias.  This is far beyond any physically feasible EDA job.
  Value = iFrag
  RunName(1:2) = 'EF'
  do Pos=8,3,-1
    Digit = int(mod(Value,36_iwp))
    RunName(Pos:Pos) = Digits(Digit+1:Digit+1)
    Value = Value/36_iwp
  end do
  if (Value /= 0) ierr = 1

end subroutine Fragment_RunName

subroutine Read_Current_RunFile(Ref,ierr,Message)

  type(EDA_Reference_Type), intent(inout) :: Ref
  integer(kind=iwp), intent(out) :: ierr
  character(len=*), intent(out) :: Message
  integer(kind=iwp) :: nSym, nData
  integer(kind=iwp) :: nBasA(1), nOrbA(1), nOccA(1)
  logical(kind=iwp) :: Found
  real(kind=wp), allocatable :: CFull(:,:)

  ierr = 1
  Message = 'Unable to read an RHF/MBPT2 reference from RunFile.'

  call EDA_Free_Reference(Ref)

  call Get_iScalar('nSym',nSym)
  if (nSym /= 1) then
    Message = 'EDA currently requires C1 (NoSym) references.'
    return
  end if

  call Get_iArray('nBas',nBasA,1)
  call Get_iArray('nOrb',nOrbA,1)
  call Get_iArray('nIsh',nOccA,1)
  call Get_iScalar('Unique atoms',Ref%nAtoms)
  call Get_cArray('Relax Method',Ref%Method,8)

  Ref%nBas = nBasA(1)
  Ref%nOrb = nOrbA(1)
  Ref%nOcc = nOccA(1)

  if ((Ref%nBas < 1) .or. (Ref%nOrb < Ref%nOcc) .or. (Ref%nOcc < 1)) then
    Message = 'Invalid basis/orbital dimensions in EDA reference RunFile.'
    return
  end if

  if ((Ref%Method(1:7) /= 'RHF-SCF') .and. (Ref%Method(1:5) /= 'MBPT2')) then
    Message = 'EDA currently accepts only RHF-SCF or closed-shell MBPT2 references.'
    return
  end if

  call qpg_dArray('SCF orbitals',Found,nData)
  if (.not. Found) then
    Message = 'SCF orbitals are missing from EDA reference RunFile.'
    return
  end if
  if (nData /= Ref%nBas*Ref%nOrb) then
    Message = 'SCF orbital array size is inconsistent with nBas*nOrb.'
    return
  end if

  allocate(CFull(Ref%nBas,Ref%nOrb))
  allocate(Ref%COcc(Ref%nBas,Ref%nOcc))
  call Get_dArray('SCF orbitals',CFull,nData)
  Ref%COcc(:,:) = CFull(:,1:Ref%nOcc)
  deallocate(CFull)

  call Get_dScalar('SCF energy',Ref%EHF)
  call Get_dScalar('Last energy',Ref%ETotal)
  Ref%HasMP2 = (Ref%Method(1:5) == 'MBPT2')
  if (.not. Ref%HasMP2) Ref%ETotal = Ref%EHF

  ierr = 0
  Message = ''

end subroutine Read_Current_RunFile

subroutine EDA_Free_Reference(Ref)

  type(EDA_Reference_Type), intent(inout) :: Ref

  if (allocated(Ref%COcc)) deallocate(Ref%COcc)
  Ref%nBas = 0
  Ref%nOrb = 0
  Ref%nOcc = 0
  Ref%nAtoms = 0
  Ref%EHF = Zero
  Ref%ETotal = Zero
  Ref%HasMP2 = .false.
  Ref%Method = '        '

end subroutine EDA_Free_Reference

end module EDA_Reference
