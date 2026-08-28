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
! Copyright (C) 2026, Yoshio Nishimoto                                 *
!***********************************************************************

#include "compiler_features.h"
#ifdef _MOLCAS_MPP_

subroutine Gather_Cholesky_Vectors(ITK,ITQ,JSYM,Array,mArray,nArray,NVLOC,NVTOT)

! Gather Cholesky vectors that are constructed in Get_Cholesky_Vectors
! The Array vector will be [NQK(ISYK=1)*NVTOT] [NQK(ISYK=2)*NVTOT]...
! NQK  : the number of (p,q) pairs per vector in the ISYK block
! NVTOT: the group's vector count over all ranks.

use CHOVEC_IO, only: NPQ_CHOTYPE
use allgather_wrapper, only: allgather
use GA_Wrapper, only: GA_NNodes, GA_NodeId
use caspt2_module, only: NSYM
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6, RtoB
use, intrinsic :: iso_c_binding, only: c_int

implicit none

integer(kind=iwp), intent(in) :: ITK, ITQ, JSYM, mArray, NVLOC, NVTOT
real(kind=wp), intent(inout) :: Array(mArray)
integer(kind=iwp), intent(out) :: nArray

integer(kind=iwp) :: ICASE, iProc, iProcEND, iProcSTA, ISYK, iVecSTA, LGLOB, LLOC, myRank, nProcs, NQK, NSEND, NVCHUNK

integer(kind=iwp), allocatable :: NVALL(:)
real(kind=wp), allocatable :: LocBuf(:)

! largest element count passing the 2 GB byte check in ALLGATHER_R
! not MAXBUF (procinp_caspt2), which bounds the ARMCI path: here MPI_Allgatherv is called directly
integer(kind=iwp), parameter :: MAXRECV = (huge(1_c_int)-mod(int(huge(1_c_int),kind=iwp),RtoB))/RtoB

! ugly hack to convert separate k/q orbital types into a specific case
ICASE = ITK*ITQ
if (ICASE == 3) then
  ICASE = 4
else
  ICASE = ICASE/2
end if

! save the local blocks (LLOC ~ mArray/nProcs elements): the larger gathered blocks overwrite them
LLOC = 0
do ISYK=1,NSYM
  LLOC = LLOC+NPQ_CHOTYPE(ICASE,ISYK,JSYM)*NVLOC
end do
call mma_allocate(LocBuf,max(LLOC,1),Label='LocBuf')
LocBuf(1:LLOC) = Array(1:LLOC)

! gather in chunks of whole processes to stay under 2 GB
! Consecutive ranks own consecutive global vector indices, so the cut points do not matter; ranks outside a chunk send nothing.
myRank = GA_NodeID()
nProcs = GA_NNodes()
call mma_allocate(NVALL,[0,nProcs-1],Label='NVALL')
NVALL(:) = 0
NVALL(myRank) = NVLOC
call GAIGOP(NVALL,nProcs,'+')
if (sum(NVALL) /= NVTOT) then
  write(u6,'(1X,A)') 'Gather_Cholesky_Vectors: local vector counts do not add up to NVTOT'
  call AbEnd()
end if

LLOC = 1
LGLOB = 1
do ISYK=1,NSYM
  NQK = NPQ_CHOTYPE(ICASE,ISYK,JSYM)
  if (NQK == 0) cycle
  iVecSTA = 0
  iProcSTA = 0
  do while (iProcSTA < nProcs)
    ! one process alone above the limit cannot be chunked any further
    if (NQK*NVALL(iProcSTA) > MAXRECV) then
      write(u6,'(1X,A,I4,A)') 'Gather_Cholesky_Vectors: the vectors of process ',iProcSTA,' exceed the 2 GB limit'
      write(u6,'(1X,A)') 'Use more processes, or a different PRHS strategy'
      call AbEnd()
    end if

    ! as many further processes as fit
    NVCHUNK = NVALL(iProcSTA)
    iProcEND = iProcSTA
    do iProc=iProcSTA+1,nProcs-1
      if (NQK*(NVCHUNK+NVALL(iProc)) > MAXRECV) exit
      NVCHUNK = NVCHUNK+NVALL(iProc)
      iProcEND = iProc
    end do
    if ((myRank >= iProcSTA) .and. (myRank <= iProcEND)) then
      NSEND = NQK*NVLOC
    else
      NSEND = 0
    end if
    call allgather(LocBuf(LLOC:),NSEND,Array(LGLOB+NQK*iVecSTA:),NQK*NVCHUNK)
    iVecSTA = iVecSTA+NVCHUNK
    iProcSTA = iProcEND+1
  end do
  LLOC = LLOC+NQK*NVLOC ! address for the next symmetry
  LGLOB = LGLOB+NQK*NVTOT
end do

nArray = LGLOB-1 ! LGLOB started at 1

call mma_deallocate(NVALL)
call mma_deallocate(LocBuf)

end subroutine Gather_Cholesky_Vectors

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(Gather_Cholesky_Vectors)

#endif
