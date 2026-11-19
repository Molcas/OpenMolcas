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

! Stored-ORDINT, C1, serial/MPI ROMP2 doubles.
! Adapted from the existing UMP2 MPI conventional development kernel.
! Real MPI: every rank reads its local distributed ORDINT contribution:
! a global sum reconstructs each AO row before cyclic row owners transform it.
! Half is then replicated by reduction.
! The second transformation/contraction is distributed over right occupied
! orbitals. With Is_Real_Par() false each displacement runs entirely locally.
! Two successive AO-pair transformations use BLAS: no Cholesky factors
! or complete four-index MO tensor are constructed. Half(q,b,j) holds
! (mu nu | j b), q=mu*(mu-1)/2+nu, mu>=nu. Each pair is unpacked to a
! symmetric AO matrix before transformation, so off-diagonal AO terms
! are included exactly twice through the two distinct matrix entries.
module ROMP2_Conventional

use Definitions, only: wp, iwp, u6, RtoB
use MBPT2_Global, only: MBPT2_All_Finite
use Constants, only: Zero, One
use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE

implicit none
private
public :: ROMP2_Conventional_Energy

contains

subroutine ROMP2_Conventional_Energy(ierr,Message)

  use Para_Info, only: nProcs, MyRank, Is_Real_Par
  use TwoDat, only: AuxTwo
  use UMP2_Global, only: nSym, nBas, nOrb, nFro, nOccA, nOccB, &
                         nOA, nOB, nVA, nVB, CAlpha, CBeta, &
                         EOccA, EOccB, EVirA, EVirB, DenTol, EAA, EBB, EAB

#include "warnings.h"

  integer(kind=iwp), intent(out) :: ierr
  character(len=*), intent(out) :: Message
  integer(kind=iwp) :: rc, CloseRC, FileNSym, FileBas(8), FileSkip(8)
  integer(kind=iwp) :: nPair, LuOrd, Workers, WorkRank
  logical(kind=iwp) :: Exists, OrdSquare, IsCho, IsDF
  logical :: Opened, Parallel, Reader, Reported
  real(kind=wp) :: Components(3)
  real(kind=wp), allocatable :: Half(:,:,:)

  ierr = 0
  Message = ''
  EAA = Zero
  EBB = Zero
  EAB = Zero
  Opened = .false.
  Reported = .false.
  LuOrd = 43
  Parallel = (nProcs > 1) .and. Is_Real_Par()
  ! ORDINT is rank-distributed in real parallel runs: every rank reads it.
  Reader = .true.
  Workers = 1
  WorkRank = 0
  if (Parallel) then
    Workers = nProcs
    WorkRank = MyRank
  end if

  call CheckReference()
  call SyncError()
  if (ierr /= 0) return
  nPair = nBas*(nBas+1)/2

  call ComputeEnergy()

  call mma_deallocate(Half,safe='*')
  if (Opened) then
    call ClsOrd(CloseRC)
    if ((CloseRC /= 0) .and. (ierr == 0)) then
      ierr = _RC_IO_ERROR_READ_
      Message = 'Unable to close ORDINT after conventional ROMP2.'
    end if
  end if
  call SyncError()
  if ((ierr == 0) .and. Parallel) then
    Components(:) = [EAA,EBB,EAB]
    call GADGOp(Components,3,'+')
    EAA = Components(1)
    EBB = Components(2)
    EAB = Components(3)
  end if
  if ((ierr == 0) .and. (.not. MBPT2_All_Finite([EAA,EBB,EAB]))) then
    ierr = _RC_GENERAL_ERROR_
    Message = 'Nonfinite globally reduced conventional ROMP2 energy.'
  end if
  if (ierr /= 0) then
    EAA = Zero
    EBB = Zero
    EAB = Zero
  end if

contains

  subroutine ComputeEnergy()
  if (Reader) then
    call OpnOrd(rc,0,'ORDINT',LuOrd)
    Opened = .true.
    if (rc /= 0) then
      ierr = _RC_IO_ERROR_READ_
      Message = 'Unable to open conventional ORDINT integrals.'
    else
      call GetOrd(rc,OrdSquare,FileNSym,FileBas,FileSkip)
      if (rc /= 0) then
        ierr = _RC_IO_ERROR_READ_
        Message = 'Unable to read the ORDINT header.'
      else if ((FileNSym /= 1) .or. (FileBas(1) /= nBas) .or. (FileSkip(1) /= 0)) then
        ierr = _RC_INPUT_ERROR_
        Message = 'ORDINT symmetry, basis dimensions, or skipped blocks disagree with the ROHF reference.'
      end if
    end if
  end if
  call SyncError()
  if (ierr /= 0) return

  ! In C1 both ordering modes expose nPair complete pair rows through
  ! RdOrd_. The ORDINT square flag governs symmetry-block duplication.

  if ((nOA >= 2) .and. (nVA >= 2)) then
    call TransformRight(CAlpha,nOccA,nOA,nVA,nOA,nVA)
    if (ierr /= 0) return
    call Contract(CAlpha,nOccA,nOA,nVA,EOccA,EVirA,nOA,nVA,EOccA,EVirA,.true.,EAA)
    call mma_deallocate(Half,safe='*')
    call SyncError()
    if (ierr /= 0) return
  end if

  if (((nOB >= 2) .and. (nVB >= 2)) .or. &
      ((nOA > 0) .and. (nVA > 0) .and. (nOB > 0) .and. (nVB > 0))) then
    ! Reuse the beta half transformation for both BB and AB.
    call TransformRight(CBeta,nOccB,nOB,nVB,max(nOA,nOB),max(nVA,nVB))
    if (ierr /= 0) return
    if ((nOB >= 2) .and. (nVB >= 2)) then
      call Contract(CBeta,nOccB,nOB,nVB,EOccB,EVirB,nOB,nVB,EOccB,EVirB,.true.,EBB)
      call SyncError()
      if (ierr /= 0) return
    end if
    if ((nOA > 0) .and. (nVA > 0)) then
      call Contract(CAlpha,nOccA,nOA,nVA,EOccA,EVirA,nOB,nVB,EOccB,EVirB,.false.,EAB)
      call SyncError()
      if (ierr /= 0) return
    end if
  end if

  end subroutine ComputeEnergy

  subroutine SyncError()
    ! Always called by every rank at matching points, even on failure.
    if (Parallel) then
      if ((ierr /= 0) .and. (.not. Reported)) &
        write(u6,'(A,I6,2A)') ' Conventional ROMP2 error on rank ',MyRank,': ',trim(Message)
      call gaIgOP_SCAL(ierr,'max')
      if (ierr /= 0) Reported = .true.
    end if
    if ((ierr /= 0) .and. (len_trim(Message) == 0)) &
      Message = 'Conventional ROMP2 failed on another MPI rank: inspect rank diagnostics.'
  end subroutine SyncError

  subroutine CheckReference()
    integer(kind=iwp) :: IntegerLimit

    if ((nSym /= 1) .or. (nBas < 1) .or. (nOrb < 1) .or. (nOrb > nBas)) then
      ierr = _RC_INPUT_ERROR_
      Message = 'Conventional ROMP2 requires a valid C1 reference.'
      return
    end if
    if ((.not. MBPT2_All_Finite(DenTol)) .or. (DenTol <= Zero)) then
      ierr = _RC_INPUT_ERROR_
      Message = 'Invalid ROMP2 denominator tolerance.'
      return
    end if
    call DecideOnCholesky(IsCho)
    call DecideOnDF(IsDF)
    if (IsCho .or. IsDF) then
      ierr = _RC_INPUT_ERROR_
      Message = 'The conventional ROMP2 reader cannot be used for Cholesky or RI integrals.'
      return
    end if
    Exists = .true.
    if (Reader) call f_Inquire('ORDINT',Exists)
    if (.not. Exists) then
      ierr = _RC_NOT_AVAILABLE_
      Message = 'Conventional ROMP2 requires stored ORDINT integrals from SEWARD: direct SCF alone is insufficient.'
      return
    end if
    if (Reader .and. AuxTwo%Opn) then
      ierr = _RC_GENERAL_ERROR_
      Message = 'Conventional ROMP2 found an already-open ordered-integral interface.'
      return
    end if
    ! Bound the product before integer multiplication, including lBuf=nPair+1.
    ! Use a variable integer limit so strict compilers do not diagnose the
    ! unavoidable loss of precision in converting HUGE(integer) to real.
    IntegerLimit = huge(nPair)
    if (real(nBas,wp)*(real(nBas,wp)+One) >= real(IntegerLimit,wp)) then
      ierr = _RC_MEMORY_ERROR_
      Message = 'Conventional ROMP2 AO-pair dimensions exceed the integer range.'
      return
    end if
  end subroutine CheckReference

  subroutine TransformRight(C,nOcc,nO,nV,MaxO,MaxV)
    integer(kind=iwp), intent(in) :: nOcc,nO,nV,MaxO,MaxV
    real(kind=wp), intent(in) :: C(nBas,nOrb)
    integer(kind=iwp) :: Available, IntegerLimit, q, b, j, nMat, iOpt, ReadRC
    real(kind=wp) :: SizeHalf, PeakFirst, PeakSecond, Required
    real(kind=wp), allocatable :: Buf(:), AO(:,:), Tmp(:,:), OV(:,:)

    ! Include both phases plus space for reader buffers and library work.
    ! Evaluate products in real arithmetic BEFORE allocating dimensions.
    SizeHalf = real(nPair,wp)*real(nV,wp)*real(nO,wp)
    PeakFirst = real(nPair,wp)+One+real(nBas,wp)**2+ &
                real(nBas,wp)*real(nO,wp)+real(nV,wp)*real(nO,wp)
    PeakSecond = real(nBas,wp)**2+real(nBas,wp)*real(MaxO,wp)+ &
                 real(MaxV,wp)*real(MaxO,wp)*real(nV,wp)
    Required = SizeHalf+max(PeakFirst,PeakSecond)
    IntegerLimit = huge(Available)
    IntegerLimit = ishft(IntegerLimit,-1)
    if (Required >= real(IntegerLimit,wp)) then
      ierr = _RC_MEMORY_ERROR_
      Message = 'Conventional ROMP2 intermediate dimensions exceed the supported integer range.'
    end if
    call SyncError()
    if (ierr /= 0) return
    call mma_maxDBLE(Available)
    if (Parallel) call gaIgOP_SCAL(Available,'min')
    if (Required+131072.0_wp > 0.8_wp*real(Available,wp)) then
      write(u6,'(A,F14.1)') ' Conventional ROMP2 estimated array storage (MiB): ',Required*real(RtoB,wp)/1048576.0_wp
      ierr = _RC_MEMORY_ERROR_
      Message = 'Insufficient memory for conventional ROMP2 half transformation: increase MOLCAS_MEM or use Cholesky.'
    end if

    call SyncError()
    if (ierr /= 0) return

    call mma_allocate(Half,nPair,nV,nO,label='ROMP2 conventional half')
    call mma_allocate(Buf,nPair+1,label='ROMP2 ORDINT row')
    call mma_allocate(AO,nBas,nBas,label='ROMP2 AO pair matrix')
    call mma_allocate(Tmp,nBas,nO,label='ROMP2 first quarter')
    call mma_allocate(OV,nV,nO,label='ROMP2 first half row')
    Half(:,:,:) = Zero
    iOpt = 1
    do q=1,nPair
      ! RdOrd_ uses (lBuf-1)/nPair to select the number of complete rows.
      ! Read the stored path directly, avoiding RdOrd's saved route flag.
      Buf(:) = Zero
      if (Reader) then
        call RdOrd_(ReadRC,iOpt,1,1,1,1,Buf,nPair+1,nMat)
        iOpt = 2
        if ((ReadRC /= 0) .or. (nMat /= 1)) then
          ierr = _RC_IO_ERROR_READ_
          Message = 'Failed to read a complete conventional AO-pair row.'
        else if (.not. MBPT2_All_Finite(Buf(1:nPair))) then
          ierr = _RC_IO_ERROR_READ_
          Message = 'Nonfinite value in stored AO integrals.'
        end if
      end if
      call SyncError()
      if (ierr /= 0) exit
      ! Sum the rank-local ORDINT contributions into the complete AO row.
      if (Parallel) call GADGOp(Buf,nPair,'+')
      if (mod(q-1,Workers) == WorkRank) then
        call Square(Buf(1:nPair),AO,nBas,1,nBas)
        call DGEMM_('N','N',nBas,nO,nBas,One,AO,nBas,C(1,nFro+1),nBas,Zero,Tmp,nBas)
        call DGEMM_('T','N',nV,nO,nBas,One,C(1,nOcc+1),nBas,Tmp,nBas,Zero,OV,nV)
        if (.not. MBPT2_All_Finite(OV)) then
          ierr = _RC_GENERAL_ERROR_
          Message = 'Nonfinite conventional half-transformed integrals.'
        end if
        do j=1,nO
          do b=1,nV
            Half(q,b,j) = OV(b,j)
          end do
        end do
      end if
      call SyncError()
      if (ierr /= 0) exit
    end do
    if ((ierr == 0) .and. Parallel) then
      ! Use contiguous pair columns, keeping each collective count bounded.
      do j=1,nO
        do b=1,nV
          call GADGOp(Half(1,b,j),nPair,'+')
        end do
      end do
    end if
    call mma_deallocate(Buf)
    call mma_deallocate(AO)
    call mma_deallocate(Tmp)
    call mma_deallocate(OV)
  end subroutine TransformRight

  subroutine Contract(C,nOcc,nO,nV,Eo,Ev,nOR,nVR,EoR,EvR,SameSpin,Energy)
    integer(kind=iwp), intent(in) :: nOcc,nO,nV,nOR,nVR
    real(kind=wp), intent(in) :: C(nBas,nOrb), Eo(nO), Ev(nV), EoR(nOR), EvR(nVR)
    logical, intent(in) :: SameSpin
    real(kind=wp), intent(out) :: Energy
    integer(kind=iwp) :: i,j,a,b
    real(kind=wp) :: Den, Numer
    real(kind=wp), allocatable :: AO(:,:), Tmp(:,:), G(:,:,:)

    Energy = Zero
    call mma_allocate(AO,nBas,nBas,label='ROMP2 second AO matrix')
    call mma_allocate(Tmp,nBas,nO,label='ROMP2 third quarter')
    call mma_allocate(G,nV,nO,nVR,label='ROMP2 occupied slice')
    ! G(a,i,b) = (i_left a_left | j_right b_right), one fixed j.
    ! Keeping all a,b for this j makes the same-spin exchange permutation
    ! G(b,i,a) available without storing a complete MO integral tensor.
    OccupiedRight: do j=1,nOR
      if (mod(j-1,Workers) /= WorkRank) cycle
      do b=1,nVR
        call Square(Half(:,b,j),AO,nBas,1,nBas)
        call DGEMM_('N','N',nBas,nO,nBas,One,AO,nBas,C(1,nFro+1),nBas,Zero,Tmp,nBas)
        call DGEMM_('T','N',nV,nO,nBas,One,C(1,nOcc+1),nBas,Tmp,nBas,Zero,G(1,1,b),nV)
      end do
      if (.not. MBPT2_All_Finite(G)) then
        ierr = _RC_GENERAL_ERROR_
        Message = 'Nonfinite conventional MO integrals.'
        exit OccupiedRight
      end if
      do i=1,nO
        if (SameSpin .and. (i == j)) cycle
        do b=1,nVR
          do a=1,nV
            if (SameSpin .and. (a == b)) cycle
            Den = Eo(i)+EoR(j)-Ev(a)-EvR(b)
            if ((.not. MBPT2_All_Finite(Den)) .or. (abs(Den) <= DenTol)) then
              ierr = _RC_GENERAL_ERROR_
              Message = 'Invalid or near-zero conventional ROMP2 orbital-energy denominator.'
              exit OccupiedRight
            end if
            if (SameSpin) then
              Numer = G(a,i,b)-G(b,i,a)
              Energy = Energy+0.25_wp*Numer*Numer/Den
            else
              Numer = G(a,i,b)
              Energy = Energy+Numer*Numer/Den
            end if
          end do
        end do
      end do
    end do OccupiedRight
    if (.not. MBPT2_All_Finite(Energy)) then
      ierr = _RC_GENERAL_ERROR_
      Message = 'Nonfinite conventional ROMP2 correlation energy.'
    end if
    call mma_deallocate(AO)
    call mma_deallocate(Tmp)
    call mma_deallocate(G)
  end subroutine Contract

end subroutine ROMP2_Conventional_Energy


end module ROMP2_Conventional
