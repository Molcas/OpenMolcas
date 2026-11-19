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

subroutine EDA(iReturn)

use Definitions, only: wp, iwp, u6
use Constants, only: auTokcalmol, Zero
use Para_Info, only: nProcs, MyRank, Is_Real_Par
use EDA_Data, only: EDA_nFrag, EDA_nAtom, EDA_RefMethod, EDA_FragBasis
use EDA_Parallel, only: EDA_Sync_Error
use EDA_Reference, only: EDA_Reference_Type, EDA_Read_Super, EDA_Read_Fragment, EDA_Free_Reference
use EDA_Energy, only: EDA_Results_Type, EDA_Compute_ClosedShell

implicit none
#include "warnings.h"

integer(kind=iwp), intent(out) :: iReturn
integer(kind=iwp) :: ierr, iFrag, nAtomsExpected, AnyError, iPass
integer :: AllocStat
character(len=256) :: Message
type(EDA_Reference_Type) :: Super
type(EDA_Reference_Type), allocatable :: Frag(:)
type(EDA_Results_Type) :: Results
logical :: IsRoot, IsParallel

IsParallel = (nProcs > 1) .and. Is_Real_Par()
IsRoot = (.not. IsParallel) .or. (MyRank == 0)

call EDA_Read_Input(ierr)
AnyError = 0
if (ierr /= _RC_ALL_IS_WELL_) AnyError = 1
if (IsParallel) call gaIgOP_SCAL(AnyError,'max')
if (AnyError /= 0) then
  iReturn = _RC_INPUT_ERROR_
  return
end if

processing: do iPass=1,1

allocate(Frag(EDA_nFrag),stat=AllocStat)
if (AllocStat /= 0) then
  Message = 'Unable to allocate EDA fragment reference data.'
  ierr = 1
  exit processing
end if

call EDA_Read_Super(Super,ierr,Message)
call EDA_Sync_Error(ierr,Message)
if (ierr /= 0) exit processing

nAtomsExpected = sum(EDA_nAtom(1:EDA_nFrag))
if (Super%nAtoms /= nAtomsExpected) then
  write(Message,'(A,I0,A,I0,A)') 'EDA NATOM sum is ',nAtomsExpected, &
                                  ' but the supermolecule RunFile contains ',Super%nAtoms,' atoms.'
  ierr = 1
else
  ierr = 0
end if
call EDA_Sync_Error(ierr,Message)
if (ierr /= 0) exit processing

if ((trim(EDA_RefMethod) == 'RHF') .and. (Super%Method /= 'RHF-SCF ')) then
  Message = 'REFERENCE=RHF requires the current supermolecule RunFile to end at RHF-SCF.'
  ierr = 1
else if ((trim(EDA_RefMethod) == 'MBPT2') .and. (.not. Super%HasMP2)) then
  Message = 'REFERENCE=MBPT2 requires the current supermolecule RunFile to contain an MBPT2 result.'
  ierr = 1
else
  ierr = 0
end if
call EDA_Sync_Error(ierr,Message)
if (ierr /= 0) exit processing

fragment_loop: do iFrag=1,EDA_nFrag
  call EDA_Read_Fragment(iFrag,Frag(iFrag),ierr,Message)
  if (ierr == 0) then
    if ((trim(EDA_RefMethod) == 'RHF') .and. (Frag(iFrag)%Method /= 'RHF-SCF ')) then
      write(Message,'(A,I0,A)') 'REFERENCE=RHF but fragment ',iFrag,' is not an RHF-SCF RunFile.'
      ierr = 1
    else if ((trim(EDA_RefMethod) == 'MBPT2') .and. (.not. Frag(iFrag)%HasMP2)) then
      write(Message,'(A,I0,A)') 'REFERENCE=MBPT2 but fragment ',iFrag,' has no MBPT2 energy.'
      ierr = 1
    end if
  end if
  call EDA_Sync_Error(ierr,Message)
  if (ierr /= 0) exit fragment_loop
end do fragment_loop
if (ierr /= 0) exit processing

call EDA_Compute_ClosedShell(Super,Frag,EDA_nFrag,EDA_FragBasis,Results,ierr,Message)
call EDA_Sync_Error(ierr,Message)
if (ierr /= 0) exit processing

if (IsRoot) call PrintResults(Super,Frag,EDA_nFrag,Results)

iReturn = _RC_ALL_IS_WELL_
end do processing

if (ierr /= 0) then
  if (IsRoot) write(u6,'(/,6X,2A)') 'EDA failed: ',trim(Message)
  iReturn = _RC_GENERAL_ERROR_
end if
call EDA_Free_Reference(Super)
if (allocated(Frag)) then
  do iFrag=1,size(Frag)
    call EDA_Free_Reference(Frag(iFrag))
  end do
  deallocate(Frag)
end if

contains

subroutine PrintResults(SuperRef,FragRef,nFrag,Res)

  type(EDA_Reference_Type), intent(in) :: SuperRef
  type(EDA_Reference_Type), intent(in) :: FragRef(:)
  integer(kind=iwp), intent(in) :: nFrag
  type(EDA_Results_Type), intent(in) :: Res
  integer(kind=iwp) :: i
  real(kind=wp) :: SumEHF

  SumEHF = Zero
  do i=1,nFrag
    SumEHF = SumEHF + FragRef(i)%EHF
  end do

  write(u6,'(/,6X,A)') '=============================================================='
  write(u6,'(6X,A)')   '  Closed-shell EDA'
  write(u6,'(6X,A)')   '=============================================================='
  write(u6,'(6X,A,1X,A)') 'Reference method:',trim(EDA_RefMethod)
  write(u6,'(6X,A,1X,A)') 'Fragment basis convention:',trim(EDA_FragBasis)
  write(u6,'(6X,A,I0)') 'Processes:',nProcs
  write(u6,'(6X,A,I0)') 'Supermolecule basis functions:',SuperRef%nBas
  write(u6,'(6X,A,I0)') 'Supermolecule occupied RHF orbitals:',SuperRef%nOcc
  write(u6,'(6X,A,F20.10,A)') 'Sum of fragment RHF energies:',SumEHF,' a.u.'
  write(u6,'(6X,A,F20.10,A)') 'Supermolecule RHF energy:   ',SuperRef%EHF,' a.u.'
  write(u6,'(/,6X,A)') 'Intermediate frozen-fragment energies'
  write(u6,'(8X,A,T54,F20.10)') 'E_X^(1): non-antisymmetric, non-orthogonal',Res%EX1
  write(u6,'(8X,A,T54,F20.10)') 'E_X^(2): antisymmetric, non-orthogonal    ',Res%EX2
  write(u6,'(8X,A,T54,F20.10)') 'E_X^(3): antisymmetric, orthogonalized    ',Res%EX3

  write(u6,'(/,6X,A)') 'EDA interaction-energy decomposition'
  write(u6,'(6X,A,T43,A,T61,A)') 'Component','Hartree','kcal/mol'
  write(u6,'(6X,A,T43,F16.8,T61,F16.4)') 'Electrostatic',Res%DeltaEele,Res%DeltaEele*auTokcalmol
  write(u6,'(6X,A,T43,F16.8,T61,F16.4)') 'Exchange',Res%DeltaEex,Res%DeltaEex*auTokcalmol
  write(u6,'(6X,A,T43,F16.8,T61,F16.4)') 'Repulsion',Res%DeltaErep,Res%DeltaErep*auTokcalmol
  write(u6,'(6X,A,T43,F16.8,T61,F16.4)') 'Polarization',Res%DeltaEpol,Res%DeltaEpol*auTokcalmol
  write(u6,'(6X,A,T43,F16.8,T61,F16.4)') 'HF interaction',Res%DeltaEHF,Res%DeltaEHF*auTokcalmol
  if (Res%HasMP2) then
    write(u6,'(6X,A,T43,F16.8,T61,F16.4)') 'MP2 dispersion',Res%DeltaEdisp,Res%DeltaEdisp*auTokcalmol
    write(u6,'(6X,A,T43,F16.8,T61,F16.4)') 'MP2 interaction',Res%DeltaEMP2,Res%DeltaEMP2*auTokcalmol
  end if
  write(u6,'(6X,A,T43,ES16.6)') 'HF closure residual',Res%Closure
  write(u6,'(6X,A/)') '=============================================================='

  ! Export the EDA quantities to the OpenMolcas verification system.
  ! A tolerance of 1.0e-7 Eh is tight enough to detect numerical
  ! regressions while allowing harmless changes in summation order.
  call Add_Info('EDA_EX1',[Res%EX1],1,7)
  call Add_Info('EDA_EX2',[Res%EX2],1,7)
  call Add_Info('EDA_EX3',[Res%EX3],1,7)
  call Add_Info('EDA_ELECTROSTATIC',[Res%DeltaEele],1,7)
  call Add_Info('EDA_EXCHANGE',[Res%DeltaEex],1,7)
  call Add_Info('EDA_REPULSION',[Res%DeltaErep],1,7)
  call Add_Info('EDA_POLARIZATION',[Res%DeltaEpol],1,7)
  call Add_Info('EDA_HF_INTERACTION',[Res%DeltaEHF],1,7)
  if (Res%HasMP2) then
    call Add_Info('EDA_MP2_DISPERSION',[Res%DeltaEdisp],1,7)
    call Add_Info('EDA_MP2_INTERACTION',[Res%DeltaEMP2],1,7)
  end if

end subroutine PrintResults

end subroutine EDA
