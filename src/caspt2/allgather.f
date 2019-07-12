************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
#ifdef _MOLCAS_MPP_

      module allgather_wrapper
      private
      public :: allgather

      interface allgather
        module procedure allgather_R, allgather_I
      end interface

      contains
      SUBROUTINE ALLGATHER_R(SEND,NSEND,RECV,NRECV)
      use mpi
      implicit none
************************************************************************
* allgather: gathers local buffers SEND of size NSEND on
*            each process into a buffer RECV of size NRECV.
*            The receiving buffer is allocated by this subroutine.
************************************************************************
#include "warnings.fh"
#include "WrkSpc.fh"

#include "global.fh"
#include "mafdecls.fh"

      real*8, intent(in) :: SEND(nSend)
      integer, intent(in) :: nSend
      real*8, intent(out) :: RECV(nRecv)
      integer, intent(in) :: nRecv

      integer*4 :: NSEND4(1), ITYPE4, IERROR4, nRecv4Tot
      INTEGER*4, ALLOCATABLE :: NRECV4(:),IDISP4(:)
      INTEGER*4, PARAMETER :: ONE4 = 1
      integer :: nBytes, myrank, nProcs, i

      ITYPE4 = MPI_REAL8
      NBYTES = 8 * NRECV

      IF (NBYTES.GT.2147483647) THEN
        WRITE(6,'(1X,A)') 'WARNING: ALLGATHER: receive buffer > 2GB'
        WRITE(6,'(1X,A)') 'some MPI implementations cannot handle this'
        WRITE(6,'(1X,A)') 'I will continue, but it might crash...'
        CALL XFLUSH(6)
      END IF

      MYRANK = GA_NODEID()
      NPROCS = GA_NNODES()

      ALLOCATE(NRECV4(0:NPROCS-1))
      ALLOCATE(IDISP4(0:NPROCS-1))

! first, gather the sendbuffer size of each process in NRECV4
      NSEND4(1)=INT(NSEND,KIND(NSEND4))
      CALL MPI_ALLGATHER(NSEND4,ONE4,MPI_INTEGER4,
     &                   NRECV4,ONE4,MPI_INTEGER4,
     &                   MPI_COMM_WORLD, IERROR4)
      IF (IERROR4.NE.0) THEN
        WRITE(6,'(1X,A,I4)') 'ERROR: ALLGATHER: MPI_Allgather ',IERROR4
        CALL ABEND()
      END IF

! check sum of send buffers against size of the receive buffer
      NRECV4TOT=0
      DO I=0,NPROCS-1
        NRECV4TOT=NRECV4TOT+NRECV4(I)
      END DO
      IF (NRECV4TOT.NE.NRECV) THEN
        WRITE(6,'(1X,A)') 'ERROR: ALLGATHER: buffer sizes do not match'
        CALL ABEND()
      END IF

! compute the displacments from the different sizes in IDISP
      IDISP4(0)=0
      DO I=1,NPROCS-1
        IDISP4(I)=IDISP4(I-1)+NRECV4(I-1)
      END DO

! gather the local send buffers into the receive buffer
      CALL MPI_ALLGATHERV(SEND,NSEND4(1),ITYPE4,
     &                    RECV,NRECV4,IDISP4,ITYPE4,
     &                    MPI_COMM_WORLD,IERROR4)
      IF (IERROR4.NE.0) THEN
        WRITE(6,'(1X,A,I4)') 'ERROR: ALLGATHER: MPI_Allgatherv ',IERROR4
        CALL ABEND()
      END IF
      end subroutine  allgather_R

      SUBROUTINE ALLGATHER_I(SEND,NSEND,RECV,NRECV)
      use mpi
      implicit none
************************************************************************
* allgather: gathers local buffers SEND of size NSEND on
*            each process into a buffer RECV of size NRECV.
*            The receiving buffer is allocated by this subroutine.
************************************************************************
#include "warnings.fh"
#include "WrkSpc.fh"

#include "global.fh"
#include "mafdecls.fh"

      integer, intent(in) :: SEND(nSend)
      integer, intent(in) :: nSend
      integer, intent(out) :: RECV(nRecv)
      integer, intent(in) :: nRecv

      integer*4 :: NSEND4(1), ITYPE4, IERROR4, nRecv4Tot
      INTEGER*4, ALLOCATABLE :: NRECV4(:),IDISP4(:)
      INTEGER*4, PARAMETER :: ONE4 = 1
      integer :: nBytes, myrank, nProcs, i

#ifdef _I8_
        ITYPE4=MPI_INTEGER8
        NBYTES=8*NRECV
#else
        ITYPE4=MPI_INTEGER4
        NBYTES=4*NRECV
#endif

      IF (NBYTES.GT.2147483647) THEN
        WRITE(6,'(1X,A)') 'WARNING: ALLGATHER: receive buffer > 2GB'
        WRITE(6,'(1X,A)') 'some MPI implementations cannot handle this'
        WRITE(6,'(1X,A)') 'I will continue, but it might crash...'
        CALL XFLUSH(6)
      END IF

      MYRANK = GA_NODEID()
      NPROCS = GA_NNODES()

      ALLOCATE(NRECV4(0:NPROCS-1))
      ALLOCATE(IDISP4(0:NPROCS-1))

! first, gather the sendbuffer size of each process in NRECV4
      NSEND4(1)=INT(NSEND,KIND(NSEND4))
      CALL MPI_ALLGATHER(NSEND4,ONE4,MPI_INTEGER4,
     &                   NRECV4,ONE4,MPI_INTEGER4,
     &                   MPI_COMM_WORLD, IERROR4)
      IF (IERROR4.NE.0) THEN
        WRITE(6,'(1X,A,I4)') 'ERROR: ALLGATHER: MPI_Allgather ',IERROR4
        CALL ABEND()
      END IF

! check sum of send buffers against size of the receive buffer
      NRECV4TOT=0
      DO I=0,NPROCS-1
        NRECV4TOT=NRECV4TOT+NRECV4(I)
      END DO
      IF (NRECV4TOT.NE.NRECV) THEN
        WRITE(6,'(1X,A)') 'ERROR: ALLGATHER: buffer sizes do not match'
        CALL ABEND()
      END IF

! compute the displacments from the different sizes in IDISP
      IDISP4(0)=0
      DO I=1,NPROCS-1
        IDISP4(I)=IDISP4(I-1)+NRECV4(I-1)
      END DO

! gather the local send buffers into the receive buffer
      CALL MPI_ALLGATHERV(SEND,NSEND4(1),ITYPE4,
     &                    RECV,NRECV4,IDISP4,ITYPE4,
     &                    MPI_COMM_WORLD,IERROR4)
      IF (IERROR4.NE.0) THEN
        WRITE(6,'(1X,A,I4)') 'ERROR: ALLGATHER: MPI_Allgatherv ',IERROR4
        CALL ABEND()
      END IF
      end subroutine allgather_I
      end module allgather_wrapper

#elif defined (NAGFOR)
c Some compilers do not like empty files
      SUBROUTINE EMPTY_ALLGATHER()
      END
#endif
