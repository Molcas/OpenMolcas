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


subroutine UMP2_Driver(ireturn)

use Definitions, only: wp, iwp, u6, RtoB
use MBPT2_Global, only: MBPT2_All_Finite
use Constants, only: Zero, Three
use Para_Info, only: nProcs, MyRank, Is_Real_Par
use Data_Structures, only: SBA_Type, Allocate_DT, Deallocate_DT
use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
use Cholesky, only: NumCho, nDimRS, ChoNSym => nSym, ChoNBas => nBas
use UMP2_Global, only: UMP2_Clean, nBas, nOccA, nOccB, nFro, &
                       nOA, nOB, nVA, nVB, EOccA, EOccB, EVirA, EVirB, &
                       DoCholesky, DenTol, ESCF, EAA, EBB, EAB, ECorr, &
                       ETotal, EnergyReady, iPL
use UMP2_Cholesky_Transform, only: UMP2_Transform_AO_Batch
use UMP2_Conventional, only: UMP2_Conventional_Energy

implicit none

interface
  subroutine RdInp(CMO,Eall,Eocc,Eext,iTst,ESCF,OpenErr,OpenMessage,SymFrozen,SymOcc)
    use Definitions, only: wp,iwp
#include "intent.fh"
    real(kind=wp),optional,intent(_OUT_) :: CMO(*),Eall(*),Eocc(*),Eext(*)
    integer(kind=iwp),optional,intent(out) :: iTst,OpenErr
    real(kind=wp),optional,intent(out) :: ESCF
    character(len=*),optional,intent(out) :: OpenMessage
    integer(kind=iwp),optional,intent(out) :: SymFrozen(:)
    integer(kind=iwp),optional,intent(in) :: SymOcc(:,:)
  end subroutine RdInp
end interface

#include "warnings.h"

integer(kind=iwp), intent(out) :: ireturn
integer(kind=iwp) :: ierr, FinalRC, nChoLocal, nChoGlobal
integer(kind=iwp) :: nRead, iVec, iRedC, GlobalVec, LocalVectorErr
integer(kind=iwp) :: ChoOffset, nChoMin, nChoMax, nWorkers, iWorker
integer(kind=iwp) :: Available, Skip(8), Dims(1), iTol, element_limit
integer(kind=iwp), external :: Cho_X_GetTol
integer(kind=iwp), allocatable :: ChoCounts(:)
logical(kind=iwp) :: IsDF, ChoReady, HaveAO
logical :: Parallel
real(kind=wp) :: NAElements, NBElements, Required
real(kind=wp), allocatable :: RedVec(:), La(:,:,:), Lb(:,:,:)
type(SBA_Type), target :: AO
character(len=256) :: Message

! MPI implementation:
!
! * Cholesky vectors remain distributed by the OpenMolcas Cholesky
!   layer. NumCho(1) is therefore the local vector count on each rank.
! * Each rank transforms only its local AO Cholesky vectors.
! * The transformed vectors are placed in disjoint rank-owned slices of
!   globally sized La/Lb arrays.
! * GADGOp reconstructs the complete transformed vector set on every
!   rank before UMP2_Energy distributes occupied-pair contractions.
!
! The global Cholesky-vector order is rank-blocked rather than the
! original serial order. This is valid because all UMP2 contractions
! depend only on complete sums over the Cholesky index P, and the same
! rank-blocked order is used for La and Lb.
!
! This routine must be entered before the original RdMBPT/RdInp.
! It exclusively owns Cholesky initialization for this MBPT2 invocation.

ireturn = _RC_GENERAL_ERROR_
ChoReady = .false.
HaveAO = .false.
Message = ''
Parallel = (nProcs > 1) .and. Is_Real_Par()

nWorkers = 1
if (Parallel) nWorkers = nProcs

call UMP2_Clean()

! Energy evaluations are also allowed at displaced geometries.
! NUMERICAL_GRADIENT disables real parallelism inside each displacement;
! the existing Is_Real_Par guards select the local serial energy path.
! Analytic UMP2 densities and gradients are not computed here.

call UMP2_Read_Reference(ierr,Message)
if (ierr /= 0) then
  call Fail(_RC_INPUT_ERROR_,Message)
  return
end if

call RdInp(OpenErr=ierr,OpenMessage=Message)
if (ierr /= 0) then
  call Fail(_RC_INPUT_ERROR_,Message)
  return
end if

call UMP2_Setup_Spaces(ierr,Message)
if (ierr /= 0) then
  call Fail(_RC_INPUT_ERROR_,Message)
  return
end if

call DecideOnCholesky(DoCholesky)
call DecideOnDF(IsDF)
if (IsDF) then
  call Fail(_RC_NOT_AVAILABLE_, &
            'UMP2 RI integrals are unsupported; use Cholesky or conventional stored AO integrals.')
  return
end if

if (.not. DoCholesky) then
  call UMP2_Conventional_Energy(ierr,Message)
  if (ierr /= 0) then
    call Fail(ierr,Message)
    return
  end if
else
  ! Initialize the OpenMolcas Cholesky interface. In a real MPI run the
  ! Cholesky layer exposes only the vectors owned by the current rank.
  call Cho_X_Init(ierr,Zero)
  if (ierr /= 0) then
    if (MyRank == 0) &
      write(u6,'(A,I8)') ' UMP2: Cho_X_Init failed, return code ',ierr
    call UMP2_Clean()
    call Quit(_RC_CHO_INI_)
    return
  end if
  ChoReady = .true.

  if ((ChoNSym /= 1) .or. (ChoNBas(1) /= nBas)) then
    call Fail(_RC_CHO_INI_, &
              'Cholesky and UHF reference dimensions disagree.')
    return
  end if

  nChoLocal = NumCho(1)
  if (nChoLocal < 0) then
    call Fail(_RC_CHO_INI_,'Invalid local Cholesky-vector count.')
    return
  end if

  ! Collect the local vector counts on every rank. The resulting rank-block
  ! layout defines the global P ordering used by both spin transformations.
  call mma_allocate(ChoCounts,nWorkers,label='UMP2 Cholesky counts')
  ChoCounts(:) = 0

  if (Parallel) then
    ChoCounts(MyRank+1) = nChoLocal
    call gaIgOP(ChoCounts(1),nWorkers,'+')
  else
    ChoCounts(1) = nChoLocal
  end if

  nChoGlobal = 0
  do iWorker=1,nWorkers
    if ((ChoCounts(iWorker) < 0) .or. &
        (ChoCounts(iWorker) > huge(nChoGlobal)-nChoGlobal)) then
      call Fail(_RC_CHO_INI_,'Global Cholesky-vector count exceeds integer range.')
      return
    end if
    nChoGlobal = nChoGlobal+ChoCounts(iWorker)
  end do

  if (nChoGlobal <= 0) then
    call Fail(_RC_CHO_INI_,'The global Cholesky-vector count is zero.')
    return
  end if

  ChoOffset = 0
  if (Parallel .and. (MyRank > 0)) then
    ChoOffset = sum(ChoCounts(1:MyRank))
  end if

  nChoMin = minval(ChoCounts)
  nChoMax = maxval(ChoCounts)

  if (.not. allocated(nDimRS)) then
    call Fail(_RC_CHO_INI_, &
              'Cholesky reduced-set dimensions are missing.')
    return
  end if

  if ((size(nDimRS,1) < 1) .or. (size(nDimRS,2) < 1)) then
    call Fail(_RC_CHO_INI_, &
              'Invalid Cholesky reduced-set dimensions.')
    return
  end if

  nRead = maxval(nDimRS(1,:))
  if (nRead <= 0) then
    call Fail(_RC_CHO_INI_,'Empty Cholesky reduced sets.')
    return
  end if

  ! La/Lb are globally replicated arrays. Check their complete
  ! per-rank memory requirement rather than the local Cholesky count.
  NAElements = real(nVA,wp)*real(nOA,wp)*real(nChoGlobal,wp)
  NBElements = real(nVB,wp)*real(nOB,wp)*real(nChoGlobal,wp)

  element_limit = huge(nRead)
  element_limit = element_limit/RtoB
  if (max(NAElements,NBElements,real(nBas,wp)**2) > real(element_limit,wp)) then
    call Fail(_RC_MEMORY_ERROR_, &
              'UMP2 array dimensions exceed the supported integer indexing range.')
    return
  end if

  Required = NAElements+NBElements+real(nRead,wp) &
             +6.0_wp*real(nBas,wp)**2 &
             +Three*real(nBas,wp)*real(max(1,nOA,nOB),wp)

  call mma_maxDBLE(Available)

  if (Required > 0.75_wp*real(Available,wp)) then
    if (MyRank == 0) then
      write(u6,'(A,ES16.6)') &
        ' Estimated UMP2 working storage (real words): ',Required
      write(u6,'(A,I16)') &
        ' Available real words:                      ',Available
    end if
    call Fail(_RC_MEMORY_ERROR_, &
              'Insufficient memory for replicated MPI UMP2 vectors.')
    return
  end if

  call mma_allocate(La,nVA,nOA,nChoGlobal,label='UMP2 global La')
  call mma_allocate(Lb,nVB,nOB,nChoGlobal,label='UMP2 global Lb')
  call mma_allocate(RedVec,nRead,label='UMP2 reduced vector')

  ! Nonlocal rank blocks must be exactly zero before the global sum.
  La(:,:,:) = Zero
  Lb(:,:,:) = Zero

  ! Case 0 gives A2(nBas*nBas,1) and A3(nBas,nBas,1) aliases.
  ! Allocate_SBA does NOT initialize ipOff: set it explicitly.
  Dims(1) = nBas
  call Allocate_DT(AO,Dims,Dims,1,1,1,0,Label='UMP2 full AO vector')
  HaveAO = .true.
  AO%ipOff = 1
  Skip(:) = 0
  Skip(1) = 1
  iRedC = -1

  if ((MyRank == 0) .and. (iPL >= 3)) then
    write(u6,'(6X,A,I10)') 'Cholesky vectors (global)              = ',nChoGlobal
    if (Parallel) then
      write(u6,'(6X,A,2I10)') &
        'Cholesky vectors/rank (min,max)       = ',nChoMin,nChoMax
    end if
  end if

  ! Each rank reads and transforms only its locally owned Cholesky vectors.
  ! The result is written directly into that rank's disjoint global slice.
  !
  ! Do not return from inside this loop in a real MPI run. A local read or
  ! transform error is reduced after the loop so all ranks reach the same
  ! collective before cleanup.
  LocalVectorErr = 0

  do iVec=1,nChoLocal

    GlobalVec = ChoOffset+iVec

    call Cho_X_getVfull(ierr,RedVec,nRead,iVec,1,1,2, &
                        iRedC,AO,Skip,.true.)

    if (ierr /= 0) then
      write(u6,'(A,I6,A,3I10)') &
        ' UMP2 rank ',MyRank,' vector/read error (local,global,rc): ', &
        iVec,GlobalVec,ierr
      LocalVectorErr = 1
      exit
    end if

    call UMP2_Transform_AO_Batch(AO%SB(1)%A3, &
                                 La(:,:,GlobalVec:GlobalVec), &
                                 Lb(:,:,GlobalVec:GlobalVec), &
                                 ierr,Message)

    if (ierr /= 0) then
      write(u6,'(A,I6,A,A)') &
        ' UMP2 rank ',MyRank,' transformation error: ',trim(Message)
      LocalVectorErr = 2
      exit
    end if

  end do

  if (Parallel) call gaIgOP_SCAL(LocalVectorErr,'max')

  if (LocalVectorErr /= 0) then
    call Fail(_RC_IO_ERROR_READ_, &
              'A local Cholesky-vector read or transformation failed.')
    return
  end if

  ! Reconstruct the complete transformed Cholesky-vector set on every rank.
  ! Each global P slice is nonzero on exactly one owning rank and zero on
  ! all other ranks, so an elementwise global sum is equivalent to an
  ! all-gather of the transformed vectors.
  if (Parallel) then

    if ((nVA > 0) .and. (nOA > 0) .and. (nChoGlobal > 0)) then
      call GADGOp(La(1,1,1),size(La),'+')
    end if

    if ((nVB > 0) .and. (nOB > 0) .and. (nChoGlobal > 0)) then
      call GADGOp(Lb(1,1,1),size(Lb),'+')
    end if

  end if

  if ((.not. MBPT2_All_Finite(La)) .or. &
      (.not. MBPT2_All_Finite(Lb))) then
    call Fail(_RC_GENERAL_ERROR_, &
              'Nonfinite globally assembled transformed Cholesky vectors.')
    return
  end if

  ! Every rank now owns the complete transformed P dimension. Therefore the
  ! Occupied-pair work sharing in UMP2_Energy is independent of the
  ! ordering of the complete Cholesky-vector set.
  call UMP2_Energy(La,Lb,EOccA,EVirA,EOccB,EVirB, &
                   nOA,nVA,nOB,nVB,nChoGlobal,DenTol, &
                   EAA,EBB,EAB,ierr)

  if (ierr /= 0) then
    if (MyRank == 0) &
      write(u6,'(A,I8)') ' UMP2 energy-kernel return code: ',ierr
    call Fail(_RC_GENERAL_ERROR_, &
              'UMP2 energy evaluation failed; check orbital-energy denominators.')
    return
  end if

end if

ECorr = EAA+EBB+EAB
ETotal = ESCF+ECorr

if (.not. MBPT2_All_Finite([ESCF,EAA,EBB,EAB,ECorr,ETotal])) then
  call Fail(_RC_GENERAL_ERROR_,'UMP2 produced a nonfinite energy.')
  return
end if

! Close files and release vector storage before publishing a result.
call ReleaseVectors(FinalRC)
if (FinalRC /= 0) then
  call Fail(_RC_CHO_RUN_,'Cholesky finalization failed.')
  return
end if

EnergyReady = .true.

if ((MyRank == 0) .and. (iPL >= 2)) then
  call PrintResultsHeader()

  write(u6,'(6X,A,T48,A)') 'Method', '= Canonical UMP2'
  if (DoCholesky) then
    write(u6,'(6X,A,T48,A)') 'Integral treatment', '= Cholesky'
  else
    write(u6,'(6X,A,T48,A)') 'Integral treatment', '= Conventional stored AO integrals'
  end if
  write(u6,'(6X,A,T48,A,I8)') 'Frozen occupied orbitals per spin','=',nFro

  if (iPL >= 3) then
    write(u6,'(6X,A,T48,A,2I8)') 'Occupied orbitals (alpha,beta)','=',nOccA,nOccB
    write(u6,'(6X,A,T48,A,2I8)') 'Correlated occupied orbitals (alpha,beta)','=',nOA,nOB
    write(u6,'(6X,A,T48,A,2I8)') 'Virtual orbitals (alpha,beta)','=',nVA,nVB
  end if

  write(u6,'(//6X,A,T48,A,F20.10,A)') 'UHF reference energy','=',ESCF,' a.u.'
  write(u6,'(6X,A,T48,A,F20.10,A)') 'Second-order correlation energy','=',ECorr,' a.u.'
  write(u6,'(6X,A,T48,A,F20.10,A)') 'Total UMP2 energy','=',ETotal,' a.u.'
  write(u6,*)
  call PrintResult(u6,'(6X,A,T50,F19.10)', &
                   'Total MBPT2 energy',0,'',[ETotal],1)
  write(u6,*)
end if

call Store_Energies(1,[ETotal],1)

! Distinguish this energy-only route from the RHF property machinery.
call Put_cArray('Relax Method','UMP2    ',8)
call Put_iScalar('mp2prpt',0)

iTol = 8
if (DoCholesky) iTol = Cho_X_GetTol(8)

call Add_Info('E_MP2',[ETotal],1,iTol)
call Add_Info('E_UMP2_AA',[EAA],1,iTol)
call Add_Info('E_UMP2_BB',[EBB],1,iTol)
call Add_Info('E_UMP2_AB',[EAB],1,iTol)

if (MyRank == 0) call xFlush(u6)

call UMP2_Clean()
ireturn = _RC_ALL_IS_WELL_

contains

subroutine PrintResultsHeader()

  write(u6,'(/,6X,A)') &
    '********************************************************************************************************'
  write(u6,'(6X,A)') &
    '*                                               Results                                                *'
  write(u6,'(6X,A/)') &
    '********************************************************************************************************'

end subroutine PrintResultsHeader


subroutine ReleaseVectors(Code)

  integer(kind=iwp), intent(out) :: Code

  Code = 0

  if (HaveAO) then
    call Deallocate_DT(AO)
    HaveAO = .false.
  end if

  call mma_deallocate(RedVec,safe='*')
  call mma_deallocate(La,safe='*')
  call mma_deallocate(Lb,safe='*')
  call mma_deallocate(ChoCounts,safe='*')

  if (ChoReady) then
    call Cho_X_Final(Code)
    ChoReady = .false.
  end if

end subroutine ReleaseVectors


subroutine Fail(Code,Text)

  integer(kind=iwp), intent(in) :: Code
  character(len=*), intent(in) :: Text

  integer(kind=iwp) :: CleanupRC

  if (MyRank == 0) &
    write(u6,'(/,A)') ' UMP2 error: '//trim(Text)

  call ReleaseVectors(CleanupRC)

  if ((CleanupRC /= 0) .and. (MyRank == 0)) &
    write(u6,'(A,I8)') ' Cholesky cleanup return code: ',CleanupRC

  call UMP2_Clean()
  ireturn = Code

end subroutine Fail

end subroutine UMP2_Driver
