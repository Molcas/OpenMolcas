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

#include "compiler_features.h"

#ifdef _MOLCAS_MPP_
module allgather_wrapper

use definitions, only: iwp, wp, u6
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
private
public :: allgather
public :: allgather_R, allgather_I

interface allgather
  module procedure :: allgather_R, allgather_I
end interface

#include "mpi_interfaces.fh"

contains

subroutine ALLGATHER_R(SEND,NSEND,RECV,NRECV)
  !*********************************************************************
  ! allgather: gathers local buffers SEND of size NSEND on
  !            each process into a buffer RECV of size NRECV.
  !            The receiving buffer is allocated by this subroutine.
  !*********************************************************************

  use mpi
  use definitions, only: MPIInt

# include "warnings.h"
# include "global.fh"
# include "mafdecls.fh"
  integer(kind=iwp), intent(in) :: nSend
  real(kind=wp), intent(in) :: SEND(nSend)
  integer(kind=iwp), intent(in) :: nRecv
  real(kind=wp), intent(out) :: RECV(nRecv)
  integer(kind=MPIInt) :: NSEND4(1), ITYPE4, IERROR4, nRecv4Tot
  integer(kind=MPIInt), allocatable :: NRECV4(:), IDISP4(:)
  integer(kind=iwp) :: nBytes, nProcs, i

  ITYPE4 = MPI_REAL8
  NBYTES = 8*NRECV

  if (NBYTES > 2147483647) then
    write(u6,'(1X,A)') 'WARNING: ALLGATHER: receive buffer > 2GB'
    write(u6,'(1X,A)') 'some MPI implementations cannot handle this'
    write(u6,'(1X,A)') 'I will continue, but it might crash...'
  end if

  NPROCS = GA_NNODES()

  call MMA_ALLOCATE(NRECV4,[0,NPROCS-1],Label='NRECV4')
  call MMA_ALLOCATE(IDISP4,[0,NPROCS-1],Label='IDISP4')

  ! first, gather the sendbuffer size of each process in NRECV4
  NSEND4(1) = int(NSEND,kind=MPIInt)
  call MPI_ALLGATHER(NSEND4,1_MPIInt,MPI_INTEGER,NRECV4,1_MPIInt,MPI_INTEGER,MPI_COMM_WORLD,IERROR4)
  if (IERROR4 /= 0) then
    write(u6,'(1X,A,I4)') 'ERROR: ALLGATHER: MPI_Allgather ',IERROR4
    call ABEND()
  end if

  ! check sum of send buffers against size of the receive buffer
  NRECV4TOT = 0
  do I=0,NPROCS-1
    NRECV4TOT = NRECV4TOT+NRECV4(I)
  end do
  if (NRECV4TOT /= NRECV) then
    write(u6,'(1X,A)') 'ERROR: ALLGATHER: buffer sizes do not match'
    call ABEND()
  end if

  ! compute the displacments from the different sizes in IDISP
  IDISP4(0) = 0
  do I=1,NPROCS-1
    IDISP4(I) = IDISP4(I-1)+NRECV4(I-1)
  end do

  ! gather the local send buffers into the receive buffer
  call MPI_ALLGATHERV(SEND,NSEND4(1),ITYPE4,RECV,NRECV4,IDISP4,ITYPE4,MPI_COMM_WORLD,IERROR4)
  if (IERROR4 /= 0) then
    write(u6,'(1X,A,I4)') 'ERROR: ALLGATHER: MPI_Allgatherv ',IERROR4
    call ABEND()
  end if
  call MMA_DEALLOCATE(NRECV4)
  call MMA_DEALLOCATE(IDISP4)

end subroutine allgather_R

subroutine ALLGATHER_I(SEND,NSEND,RECV,NRECV)
  !*********************************************************************
  ! allgather: gathers local buffers SEND of size NSEND on
  !            each process into a buffer RECV of size NRECV.
  !            The receiving buffer is allocated by this subroutine.
  !*********************************************************************

  use mpi
  use definitions, only: MPIInt

#include "warnings.h"
#include "global.fh"
#include "mafdecls.fh"
  integer(kind=iwp), intent(in) :: nSend
  integer(kind=iwp), intent(in) :: SEND(nSend)
  integer(kind=iwp), intent(in) :: nRecv
  integer(kind=iwp), intent(out) :: RECV(nRecv)
  integer(kind=MPIInt) :: NSEND4(1), ITYPE4, IERROR4, nRecv4Tot
  integer(kind=MPIInt), allocatable :: NRECV4(:), IDISP4(:)
  integer(kind=iwp) :: nBytes, nProcs, i

# ifdef _I8_
  ITYPE4 = MPI_INTEGER8
  NBYTES = 8*NRECV
# else
  ITYPE4 = MPI_INTEGER4
  NBYTES = 4*NRECV
# endif

  if (NBYTES > 2147483647) then
    write(u6,'(1X,A)') 'WARNING: ALLGATHER: receive buffer > 2GB'
    write(u6,'(1X,A)') 'some MPI implementations cannot handle this'
    write(u6,'(1X,A)') 'I will continue, but it might crash...'
  end if

  NPROCS = GA_NNODES()

  call MMA_ALLOCATE(NRECV4,[0,NPROCS-1],Label='NRECV4')
  call MMA_ALLOCATE(IDISP4,[0,NPROCS-1],Label='IDISP4')

  ! first, gather the sendbuffer size of each process in NRECV4
  NSEND4(1) = int(NSEND,kind=MPIInt)
  call MPI_ALLGATHER(NSEND4,1_MPIInt,MPI_INTEGER,NRECV4,1_MPIInt,MPI_INTEGER,MPI_COMM_WORLD,IERROR4)
  if (IERROR4 /= 0) then
    write(u6,'(1X,A,I4)') 'ERROR: ALLGATHER: MPI_Allgather ',IERROR4
    call ABEND()
  end if

  ! check sum of send buffers against size of the receive buffer
  NRECV4TOT = 0
  do I=0,NPROCS-1
    NRECV4TOT = NRECV4TOT+NRECV4(I)
  end do
  if (NRECV4TOT /= NRECV) then
    write(u6,'(1X,A)') 'ERROR: ALLGATHER: buffer sizes do not match'
    call ABEND()
  end if

  ! compute the displacments from the different sizes in IDISP
  IDISP4(0) = 0
  do I=1,NPROCS-1
    IDISP4(I) = IDISP4(I-1)+NRECV4(I-1)
  end do

  ! gather the local send buffers into the receive buffer
  call MPI_ALLGATHERV(SEND,NSEND4(1),ITYPE4,RECV,NRECV4,IDISP4,ITYPE4,MPI_COMM_WORLD,IERROR4)
  if (IERROR4 /= 0) then
    write(u6,'(1X,A,I4)') 'ERROR: ALLGATHER: MPI_Allgatherv ',IERROR4
    call ABEND()
  end if
  call MMA_DEALLOCATE(NRECV4)
  call MMA_DEALLOCATE(IDISP4)
end subroutine allgather_I

end module allgather_wrapper

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
subroutine empty_ALLGATHER()
end subroutine empty_ALLGATHER

#endif
