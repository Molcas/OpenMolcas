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

subroutine GATHER_STRIPED(lg_M,NROW,NCOL,JSTA,NCOLB,BLK)

! Gather the columns JSTA to JSTA+NCOLB-1 of a horizontally striped global array into BLK(NROW,NCOLB), replicated on every process
! The buffer is then NROW*NCOLB per process, not NROW*NCOL (the full matrix)
!
! Motivation: calling GA_DGEMM directly is too slow, so replace with local DGEMM

use GA_Wrapper, only: DBL_MB, GA_NNodes, GA_NodeId
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: lg_M, NROW, NCOL, JSTA, NCOLB
real(kind=wp), intent(out) :: BLK(NROW,NCOLB)

integer(kind=iwp) :: iHi, iLo, iOff, iProc, J, jHi, jLo, LDM, mM, myRank, nLoc, nProcs, nRowP
real(kind=wp) :: DUMMY(1)

real(kind=wp), allocatable :: RECV(:)

if ((JSTA < 1) .or. ((JSTA+NCOLB-1) > NCOL)) then
  write(u6,'(1X,A)') 'GATHER_STRIPED: column range outside the array, ABORT'
  call Abend()
end if

myRank = GA_NodeID()
nProcs = GA_NNodes()

call GA_Distribution(lg_M,myRank,iLo,iHi,jLo,jHi)
nLoc = 0
if (iLo /= 0) nLoc = (iHi-iLo+1)*NCOLB

call mma_allocate(RECV,NROW*NCOLB,Label='RECV')

if (iLo /= 0) then
  if ((jHi-jLo+1) /= NCOL) then
    write(u6,'(1X,A)') 'GATHER_STRIPED: array is not striped by rows, ABORT'
    call Abend()
  end if
  call GA_Access(lg_M,iLo,iHi,jLo,jHi,mM,LDM)
  ! a row-striped block holds all columns and is stored contiguously; the same assumption is made by PSBMAT_WRITE
  if (LDM /= (iHi-iLo+1)) then
    write(u6,'(1X,A)') 'GATHER_STRIPED: unexpected leading dimension, ABORT'
    call Abend()
  end if
  ! The local patch is passed by argument association on purpose: DBL_MB is declared with size 2 (mafdecls.fh),
  ! so a section of it is out of bounds and cannot be handed to the generic ALLGATHER interface (see RHS_DAXPY).
  call GATHER_STRIPED_SEND(DBL_MB(mM+(JSTA-1)*LDM),nLoc,RECV,NROW*NCOLB)
  call GA_Release(lg_M,iLo,iHi,jLo,jHi)
else
  DUMMY(1) = Zero
  call GATHER_STRIPED_SEND(DUMMY,nLoc,RECV,NROW*NCOLB)
end if

! The gathered buffer holds one (nRowP x NCOLB) block per process in rank order,
! so it has to be scattered into the column-major result.
iOff = 0
do iProc=0,nProcs-1
  call GA_Distribution(lg_M,iProc,iLo,iHi,jLo,jHi)
  if (iLo == 0) cycle
  nRowP = iHi-iLo+1
  do J=1,NCOLB
    BLK(iLo:iHi,J) = RECV(iOff+(J-1)*nRowP+1:iOff+J*nRowP)
  end do
  iOff = iOff+nRowP*NCOLB
end do

call mma_deallocate(RECV)

contains

!-----------------------------------------------------------------------

subroutine GATHER_STRIPED_SEND(SEND,NSEND,RECV,NRECV)

use allgather_wrapper, only: allgather
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NSEND, NRECV
real(kind=wp), intent(in) :: SEND(NSEND)
real(kind=wp), intent(out) :: RECV(NRECV)

call allgather(SEND,NSEND,RECV,NRECV)

end subroutine GATHER_STRIPED_SEND

end subroutine GATHER_STRIPED

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(GATHER_STRIPED)

#endif
