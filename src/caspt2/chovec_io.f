************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Steven Vancoillie                                      *
************************************************************************
      MODULE CHOVEC_IO
      use definitions, only: iwp, wp
C SVC: subroutines to read/write transformed cholesky vectors from/to
C disk. These are used in tracho (where they are written) and in rhsod
C (where they are read into a global array).
C
C Cholesky vectors are distributed in batches and transformed as
C such. They are written the same way to disk:
C -----------------------------------------------------
C |IB=1               |IB=2 |... |      |IB=NBATCH_TOT|
C -----------------------------------------------------
C |Cholesky vectors   |     |    |      |   ...       |
C |(p,q,nvec_chobatch)|     |    |      |   ...       |
C |for all itk>itq    |     |    |      |   ...       |
C |and isyk,isyq=jsym |     |    |      |   ...       |
C |for jsym=1         |     |    |      | jsym=NSYM   |
C -----------------------------------------------------
C each batch contains blocks of a certain type (case and symmetry).
C Always all vectors are present, with the elements distributed over
C batches as well as processes.
C
C Cases are enumerated as follows:
C   1 =    active,inactive
C   2 =    active,active
C   3 = secondary,active
C   4 = secondary,inactive
C
C Symmetry is determined by total symmetry JSYM and the symmetry of one
C of the orbital cases ISYQ, thus all values of Mul(ISYQ,JSYM).
C
C When reading the vectors, we group all batches of the same type
C together, so that we have all vectors of one type as one block,
C as this is how they are used to compute the integrals for RHS.

      ! indexing of distributed cholesky vectors
      !-----------------------------------------
      ! The cholesky vectors are distributed over batches and processed
      ! like that by tracho2. A batch has a specified number of cholesky
      ! vectors (batches for different JSYM are different, so are the
      ! number of vectors in them). The offset of a batch on disk is
      ! then determined by the type (case,isyq,jsym) and total batch
      ! number (ib).
      INTEGER, ALLOCATABLE, SAVE :: NVLOC_CHOBATCH(:)
      INTEGER, ALLOCATABLE, SAVE :: IDLOC_CHOGROUP(:,:,:,:)

      ! indexing of collected cholesky vectors
      !---------------------------------------
      ! The cholesky vectors can be collected together on disk from
      ! different processes. The combined sizes are stored in separate
      ! arrays.
      INTEGER(KIND=IWP), ALLOCATABLE, SAVE :: NVGLB_CHOBATCH(:)
      INTEGER(KIND=IWP), ALLOCATABLE, SAVE :: IDGLB_CHOGROUP(:,:,:,:)

      ! total amount of cholesky vectors in a certain symmetry
      INTEGER(KIND=IWP), SAVE :: NVTOT_CHOSYM(8)

      CONTAINS

************************************************************************
      FUNCTION NPQ_CHOTYPE(ICASE,ISYQ,JSYM)
************************************************************************
* Compute the number of orbital pairs for a given case (valid pair of
* inactive,active,secondary), total symmetry JSYM and component symmetry
* ISYQ. Note that the component Q which determines the symmetry block is
* the _lower_ orbital partition (e.g. inactive for active,inactive),
* which is also the slowest varying index of the pair P,Q.
************************************************************************
      use Symmetry_Info, only: Mul
      use caspt2_module, only: nAsh, nIsh, nSSh
      IMPLICIT NONE

      INTEGER(KIND=IWP) NPQ_CHOTYPE

      INTEGER(KIND=IWP), INTENT(IN) :: ICASE,ISYQ,JSYM

      INTEGER(KIND=IWP) :: ISYP,NP,NQ

      ISYP=Mul(ISYQ,JSYM)
      SELECT CASE(ICASE)
      CASE(1)
        NP=NASH(ISYP)
        NQ=NISH(ISYQ)
      CASE(2)
        NP=NASH(ISYP)
        NQ=NASH(ISYQ)
      CASE(3)
        NP=NSSH(ISYP)
        NQ=NASH(ISYQ)
      CASE(4)
        NP=NSSH(ISYP)
        NQ=NISH(ISYQ)
      CASE DEFAULT
        NP=0
        NQ=0
        CALL SYSABENDMSG('NPQ_CHOTYPE',
     &    'invalid case number', '')
      END SELECT
      NPQ_CHOTYPE=NP*NQ
      END FUNCTION NPQ_CHOTYPE

************************************************************************
      SUBROUTINE CHOVEC_SIZE(ICASE,NCHOBUF,IOFF)
************************************************************************
* Allocate a buffer to hold all cholesky vectors of type ITK,ITQ
************************************************************************
      use Symmetry_Info, only: Mul
      use caspt2_module, only: nSym
      IMPLICIT NONE

      INTEGER(KIND=IWP), INTENT(IN) :: ICASE
      INTEGER(KIND=IWP), INTENT(OUT) :: NCHOBUF,IOFF(8,8)

      INTEGER(KIND=IWP) :: ISYK,ISYQ,JSYM
      INTEGER(KIND=IWP) :: NPQ,NVTOT

      NCHOBUF=0
      DO JSYM=1,NSYM
        NVTOT=NVTOT_CHOSYM(JSYM)
        DO ISYQ=1,NSYM
          ISYK=Mul(ISYQ,JSYM)
          IOFF(ISYK,ISYQ)=NCHOBUF
          NPQ=NPQ_CHOTYPE(ICASE,ISYQ,JSYM)
          NCHOBUF=NCHOBUF+NPQ*NVTOT
        END DO
      END DO

      END SUBROUTINE CHOVEC_SIZE

************************************************************************
      SUBROUTINE CHOVEC_READ(ICASE,CHOBUF,nCHOBUF)
************************************************************************
* Read (transposed) cholesky vectors from disk, they
* are indexed as CHOBUF(IVEC,IQ,IK)
************************************************************************
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
      use caspt2_global, only: LUDRATOT
#endif
      use caspt2_global, only: LUDRA
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: nSym, nBtches, nBtch
      IMPLICIT NONE
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif
      INTEGER(KIND=IWP), INTENT(IN):: ICASE, nCHOBUF
      REAL*8, INTENT(INOUT):: CHOBUF(nCHOBUF)

      INTEGER(KIND=IWP) :: I,J,IOFF,IDISK
      INTEGER(KIND=IWP) :: IB,IBSTA,IBEND,IBOFF
      INTEGER(KIND=IWP) :: JSYM,ISYQ
      INTEGER(KIND=IWP) :: NBUF,NPQ,NV,NVTOT
      REAL*8, ALLOCATABLE:: BUF(:)

      IOFF=0
      DO JSYM=1,NSYM
        NVTOT=NVTOT_CHOSYM(JSYM)
        IBSTA=NBTCHES(JSYM)+1
        IBEND=NBTCHES(JSYM)+NBTCH(JSYM)
        DO ISYQ=1,NSYM
          NPQ=NPQ_CHOTYPE(ICASE,ISYQ,JSYM)
          IBOFF=0
          DO IB=IBSTA,IBEND
            NV=NVGLB_CHOBATCH(IB)
            NBUF=NPQ*NV
            CALL mma_allocate(BUF,NBUF,LABEL='BUF')
            IDISK=IDGLB_CHOGROUP(ICASE,ISYQ,JSYM,IB)
#ifdef _MOLCAS_MPP_
            IF (Is_Real_Par()) THEN
              ! cholesky vectors already transposed
              CALL DDAFILE(LUDRATOT,2,BUF,NBUF,IDISK)
              DO J=1,NPQ
                DO I=1,NV
                  CHOBUF(IOFF+IBOFF+I+NVTOT*(J-1))=
     &            BUF(I+NV*(J-1))
                END DO
              END DO
            ELSE
#endif
              ! cholesky vectors not transposed
              CALL DDAFILE(LUDRA,2,BUF,NBUF,IDISK)
              DO J=1,NPQ
                DO I=1,NV
                  CHOBUF(IOFF+IBOFF+I+NVTOT*(J-1))=
     &            BUF(J+NPQ*(I-1))
                END DO
              END DO
#ifdef _MOLCAS_MPP_
            ENDIF
#endif
            CALL mma_deallocate(BUF)
            IBOFF=IBOFF+NV
          END DO
          IOFF=IOFF+NVTOT*NPQ
        END DO
      END DO

      END SUBROUTINE CHOVEC_READ

************************************************************************
      SUBROUTINE CHOVEC_SAVE(CHOBUF,ICASE,ISYQ,JSYM,IB)
************************************************************************
* Write Cholesky vectors to disk.
************************************************************************
      use caspt2_global, only: LUDRA
#ifdef _DEBUGPRINT_
      use definitions, only: u6
#endif
      IMPLICIT NONE
#include "warnings.h"
      INTEGER(KIND=IWP), INTENT(IN):: ICASE,ISYQ,JSYM,IB
      REAL(KIND=WP), INTENT(INOUT):: CHOBUF(*)

      INTEGER(KIND=IWP) NPQ, JNUM, IDISK
#ifdef _DEBUGPRINT_
      INTEGER(KIND=IWP) NBUF
      REAL(KIND=WP) SQFP
      REAL(KIND=WP), EXTERNAL:: DNRM2_
#endif

C always write the chunks to LUDRA, both for serial and parallel
      NPQ=NPQ_CHOTYPE(ICASE,ISYQ,JSYM)
      JNUM=NVLOC_CHOBATCH(IB)
      IDISK=IDLOC_CHOGROUP(ICASE,ISYQ,JSYM,IB)
      CALL DDAFILE(LUDRA,1,CHOBUF,NPQ*JNUM,IDISK)

#ifdef _DEBUGPRINT_
      NBUF=NPQ*JNUM
      SQFP = DNRM2_(NBUF,CHOBUF,1)
      WRITE(u6,'(1X,A,I9,A,A,I2,A,A,I2,A,A,I2,A,A,F21.14)')
     &  'BATCH ',IB,   ', ',
     &  'CASE ' ,ICASE,', ',
     &  'ISYQ ' ,ISYQ, ', ',
     &  'JSYM ' ,JSYM, ', ',
     &  'DNRM2 ',SQFP
#endif
      END SUBROUTINE CHOVEC_SAVE

************************************************************************
      SUBROUTINE CHOVEC_LOAD(CHOBUF,ICASE,ISYQ,JSYM,IB)
************************************************************************
* Read Cholesky vectors from disk.
************************************************************************
      use caspt2_global, only: LUDRA
#ifdef _DEBUGPRINT_
      use definitions, only: u6
#endif
      IMPLICIT NONE
#include "warnings.h"
      INTEGER(KIND=IWP), INTENT(IN):: ICASE,ISYQ,JSYM,IB
      REAL(KIND=WP), INTENT(OUT):: CHOBUF(*)

      INTEGER(KIND=IWP) NPQ, JNUM, IDISK
#ifdef _DEBUGPRINT_
      INTEGER(KIND=IWP) NBUF
      REAL(KIND=WP) SQFP
      REAL(KIND=WP), EXTERNAL:: DNRM2_
#endif
C always write the chunks to LUDRA, both for serial and parallel
      NPQ=NPQ_CHOTYPE(ICASE,ISYQ,JSYM)
      JNUM=NVLOC_CHOBATCH(IB)
      IDISK=IDLOC_CHOGROUP(ICASE,ISYQ,JSYM,IB)
      CALL DDAFILE(LUDRA,2,CHOBUF,NPQ*JNUM,IDISK)

#ifdef _DEBUGPRINT_
      NBUF=NPQ*JNUM
      SQFP = DNRM2_(NBUF,CHOBUF,1)
      WRITE(u6,'(1X,A,I9,A,A,I2,A,A,I2,A,A,I2,A,A,F21.14)')
     &  'BATCH ',IB,   ', ',
     &  'CASE ' ,ICASE,', ',
     &  'ISYQ ' ,ISYQ, ', ',
     &  'JSYM ' ,JSYM, ', ',
     &  'DNRM2 ',SQFP
#endif
      END SUBROUTINE CHOVEC_LOAD

************************************************************************
      SUBROUTINE CHOVEC_COLL(CHOBUF,ICASE,ISYQ,JSYM,IB)
************************************************************************
* Routine to gather locally available cholesky vectors and collect
* all of them on each process in case of parallel run.
************************************************************************
#ifdef _MOLCAS_MPP_
      USE MPI, only: MPI_REAL8, MPI_COMM_WORLD
#  ifdef _I8_
      USE MPI, only: MPI_INTEGER8, MPI_INTEGER
#  else
      USE MPI, only: MPI_INTEGER4, MPI_INTEGER
#endif
      USE Para_Info, ONLY: nProcs, Is_Real_Par
      use caspt2_global, only: LUDRATOT
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: MPIInt
      use chocaspt2, only: NFTSPC_TOT
      use caspt2_module, only: RHSDirect
#endif
      IMPLICIT NONE
#include "warnings.h"
      REAL(KIND=WP), INTENT(INOUT) :: CHOBUF(*)
      INTEGER(KIND=IWP), INTENT(IN) :: ICASE,ISYQ,JSYM,IB

#ifdef _MOLCAS_MPP_
#  include "global.fh"
#  include "mafdecls.fh"
      integer(kind=MPIInt) IERROR4,ITYPE
      integer(kind=MPIInt), PARAMETER :: ONE4 = 1
      INTEGER :: I,JNUM,JNUMT,NPQ,NUMSEND(1),IDISKT,IERROR
      INTEGER(kind=MPIInt), ALLOCATABLE:: DISP(:), SIZE(:)
      REAL(KIND=WP), ALLOCATABLE:: TRANSP(:), RECVBUF(:)
#ifdef _DEBUGPRINT_
      INTEGER(KIND=IWP) :: MY_N,NOFF
      REAL(KIND=WP) :: SQFP
      REAL(KIND=WP), EXTERNAL :: DDOT_
#endif
#endif

#ifdef _MOLCAS_MPP_
#  ifdef _I8_
      ITYPE=MPI_INTEGER8
#  else
      ITYPE=MPI_INTEGER4
#  endif
      ITYPE=MPI_INTEGER
      IF (Is_Real_Par()) THEN
C for true parallel, also communicate chunks to each process, write them
C to LUDRATOT, so first allocate memory for the fully transformed
C vectors, and for the per-process size and offset into LUDRATOT
        CALL mma_allocate(DISP,NPROCS,Label='DISP')
        CALL mma_allocate(SIZE,NPROCS,LABEL='SIZE')

C gather sizes of local cholesky bits
        NPQ=NPQ_CHOTYPE(ICASE,ISYQ,JSYM)
        JNUM=NVLOC_CHOBATCH(IB)
        NUMSEND(1)=NPQ*JNUM
        CALL MPI_Allgather(NUMSEND,ONE4,ITYPE,SIZE(1:NPROCS),ONE4,ITYPE,
     &                     MPI_COMM_WORLD, IERROR4)
C compute offsets into the receiving array
        DISP(1)=0
        DO I=2,NPROCS
          DISP(I)=DISP(I-1)+SIZE(I-1)
        END DO

C collect the vectors
        CALL mma_allocate(RECVBUF,NFTSPC_TOT,Label='RECVBUF')
        CALL MPI_Barrier(MPI_COMM_WORLD, IERROR4)
        CALL MPI_Allgatherv_(CHOBUF,NUMSEND(1),MPI_REAL8,
     &                       RECVBUF,SIZE,DISP,
     &                       MPI_REAL8,MPI_COMM_WORLD, IERROR)

        JNUMT=NVGLB_CHOBATCH(IB)
        ! disk offset is block offset + preceding block size
        IDISKT=IDGLB_CHOGROUP(ICASE,ISYQ,JSYM,IB)

CSVC: for RHS on demand, write transposed chovecs, else just write
        IF (RHSDIRECT) THEN
          CALL mma_allocate(TRANSP,NPQ*JNUMT,Label='TRANSP')
          CALL DTRANS(NPQ,JNUMT,RECVBUF,NPQ,TRANSP,JNUMT)
          CALL DDAFILE(LUDRATOT,1,TRANSP,NPQ*JNUMT,IDISKT)
          CALL mma_deallocate(TRANSP)
        ELSE
          CALL DDAFILE(LUDRATOT,1,RECVBUF,NPQ*JNUMT,IDISKT)
        END IF

#  ifdef _DEBUGPRINT_
        WRITE(u6,*) ' process block, size, offset, fingerprint'
        DO I=1,NPROCS
          MY_N = SIZE(I)
          NOFF = 1+DISP(I)
          SQFP =DDOT_(MY_N,RECVBUF(NOFF:),1,RECVBUF(NOFF:),1)
          WRITE(u6,'(A,I6,A,2I12,ES20.12)') ' [',I,'] ',MY_N,NOFF,SQFP
        END DO
#  endif
        CALL mma_deallocate(RECVBUF)
        CALL mma_deallocate(DISP)
        CALL mma_deallocate(SIZE)
      END IF
#else
C Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_real_array(CHOBUF)
        CALL Unused_integer(ICASE)
        CALL Unused_integer(ISYQ)
        CALL Unused_integer(JSYM)
        CALL Unused_integer(IB)
      END IF
#endif
      END SUBROUTINE CHOVEC_COLL

#ifdef _MOLCAS_MPP_
************************************************************************
      SUBROUTINE MPI_Allgatherv_(SENDBUF,NSEND,MPITYPES,
     &                     RCVBUF,NRCV,NOFF,MPITYPER,MPICOMM,IERROR)
************************************************************************
* Wrapper to MPI_Allgatherv dealing with ILP64 incompatibility.
************************************************************************
      USE MPI, only: MPI_COMM_WORLD
      use definitions, only: MPIInt
      IMPLICIT NONE
      REAL(KIND=WP), INTENT(INOUT):: SENDBUF(*)
      INTEGER(KIND=IWP), INTENT(IN):: NSEND
      integer(kind=MPIInt), INTENT(IN) :: MPITYPES
      REAL(KIND=WP), INTENT(INOUT):: RCVBUF(*)
      integer(kind=MPIInt), INTENT(IN) :: NRCV(*), NOFF(*)
      integer(kind=MPIInt), INTENT(IN) :: MPITYPER, MPICOMM
      INTEGER(KIND=IWP), INTENT(OUT) :: IERROR

      integer(kind=MPIInt) :: NPROCS
      integer(kind=MPIInt) :: NSEND4
      integer(kind=MPIInt),ALLOCATABLE :: NRCV4(:),NOFF4(:)
      integer(kind=MPIInt) :: IERROR4
      INTEGER(KIND=IWP), PARAMETER :: I4=KIND(NSEND4)

      INTEGER(KIND=IWP) :: I
#ifdef _I8_
      INTEGER(KIND=IWP) :: NRCVTOT
#endif

      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS,IERROR4)

#ifdef _I8_
      NRCVTOT=0
      DO I=1,NPROCS
        NRCVTOT=NRCVTOT+NRCV(I)
      END DO
      IF (8*NRCVTOT.GT.2147483647) THEN
        WRITE(6,'(1X,A)') 'MPI_Allgatherv: total rcv buf > 2**31-1'
        WRITE(6,'(1X,A)') 'workaround to avoid buffers >2GB failed'
        WRITE(6,'(1X,A)') '-> please report this as a bug!'
        CALL ABEND()
      END IF
#endif

      ALLOCATE(NRCV4(NPROCS))
      ALLOCATE(NOFF4(NPROCS))
      NSEND4=INT(NSEND,I4)
      DO I=1,NPROCS
        NRCV4(I)=INT(NRCV(I),I4)
        NOFF4(I)=INT(NOFF(I),I4)
      END DO
      CALL MPI_Allgatherv(SENDBUF,NSEND4,MPITYPES,
     &                    RCVBUF,NRCV4,NOFF4,MPITYPER,
     &                    MPICOMM,IERROR4)

      IERROR=IERROR4
      END SUBROUTINE MPI_Allgatherv_
#endif

      END MODULE CHOVEC_IO
