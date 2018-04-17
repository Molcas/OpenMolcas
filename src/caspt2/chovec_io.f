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
C of the orbital cases ISYQ, thus all values of MUL(ISYQ,JSYM).
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
      INTEGER, ALLOCATABLE, SAVE :: NVGLB_CHOBATCH(:)
      INTEGER, ALLOCATABLE, SAVE :: IDGLB_CHOGROUP(:,:,:,:)

      ! total amount of cholesky vectors in a certain symmetry
      INTEGER, SAVE :: NVTOT_CHOSYM(8)

      CONTAINS

************************************************************************
      INTEGER FUNCTION NPQ_CHOTYPE(ICASE,ISYQ,JSYM)
************************************************************************
* Compute the number of orbital pairs for a given case (valid pair of
* inactive,active,secondary), total symmetry JSYM and component symmetry
* ISYQ. Note that the component Q which determines the symmetry block is
* the _lower_ orbital partition (e.g. inactive for active,inactive),
* which is also the slowest varying index of the pair P,Q.
************************************************************************
      IMPLICIT NONE
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"

      INTEGER :: ICASE,ISYQ,JSYM
      INTEGER :: ISYP,NP,NQ

      ISYP=MUL(ISYQ,JSYM)
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
      END FUNCTION

************************************************************************
      SUBROUTINE CHOVEC_SIZE(ICASE,NCHOBUF,IOFF)
************************************************************************
* Allocate a buffer to hold all cholesky vectors of type ITK,ITQ
************************************************************************
      IMPLICIT NONE
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"

      INTEGER :: ICASE,NCHOBUF,IOFF(8,8)
      INTEGER :: ISYK,ISYQ,JSYM
      INTEGER :: NPQ,NVTOT

      NCHOBUF=0
      DO JSYM=1,NSYM
        NVTOT=NVTOT_CHOSYM(JSYM)
        DO ISYQ=1,NSYM
          ISYK=MUL(ISYQ,JSYM)
          IOFF(ISYK,ISYQ)=NCHOBUF
          NPQ=NPQ_CHOTYPE(ICASE,ISYQ,JSYM)
          NCHOBUF=NCHOBUF+NPQ*NVTOT
        END DO
      END DO

      END SUBROUTINE

************************************************************************
      SUBROUTINE CHOVEC_READ(ICASE,LCHOBUF)
************************************************************************
* Read (transposed) cholesky vectors from disk, they
* are indexed as CHOBUF(IVEC,IQ,IK)
************************************************************************
      IMPLICIT NONE
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
      LOGICAL :: IS_REAL_PAR
#endif

      INTEGER :: ICASE,LCHOBUF

      INTEGER :: I,J,IOFF,IDISK
      INTEGER :: IB,IBSTA,IBEND,IBOFF
      INTEGER :: JSYM,ISYQ
      INTEGER :: LBUF,NBUF,NPQ,NV,NVTOT

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
            CALL GETMEM('BUF','ALLO','REAL',LBUF,NBUF)
            IDISK=IDGLB_CHOGROUP(ICASE,ISYQ,JSYM,IB)
#ifdef _MOLCAS_MPP_
            IF (Is_Real_Par()) THEN
              ! cholesky vectors already transposed
              CALL DDAFILE(LUDRATOT,2,WORK(LBUF),NBUF,IDISK)
              DO J=1,NPQ
                DO I=1,NV
                  WORK(LCHOBUF+IOFF+IBOFF+I-1+NVTOT*(J-1))=
     &            WORK(LBUF+I-1+NV*(J-1))
                END DO
              END DO
            ELSE
              ! cholesky vectors not transposed
              CALL DDAFILE(LUDRA,2,WORK(LBUF),NBUF,IDISK)
              DO J=1,NPQ
                DO I=1,NV
                  WORK(LCHOBUF+IOFF+IBOFF+I-1+NVTOT*(J-1))=
     &            WORK(LBUF+J-1+NPQ*(I-1))
                END DO
              END DO
            ENDIF
#else
            ! cholesky vectors not transposed
            CALL DDAFILE(LUDRA,2,WORK(LBUF),NBUF,IDISK)
            DO J=1,NPQ
              DO I=1,NV
                WORK(LCHOBUF+IOFF+IBOFF+I-1+NVTOT*(J-1))=
     &            WORK(LBUF+J-1+NPQ*(I-1))
              END DO
            END DO
#endif
            CALL GETMEM('BUF','FREE','REAL',LBUF,NBUF)
            IBOFF=IBOFF+NV
          END DO
          IOFF=IOFF+NVTOT*NPQ
        END DO
      END DO

      END SUBROUTINE

************************************************************************
      SUBROUTINE CHOVEC_SAVE(CHOBUF,ICASE,ISYQ,JSYM,IB)
************************************************************************
* Write Cholesky vectors to disk.
************************************************************************
      Implicit real*8 (a-h,o-z)
#include "rasdim.fh"
#include "warnings.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "chocaspt2.fh"
      DIMENSION CHOBUF(*)

C always write the chunks to LUDRA, both for serial and parallel
      NPQ=NPQ_CHOTYPE(ICASE,ISYQ,JSYM)
      JNUM=NVLOC_CHOBATCH(IB)
      IDISK=IDLOC_CHOGROUP(ICASE,ISYQ,JSYM,IB)
      CALL DDAFILE(LUDRA,1,CHOBUF,NPQ*JNUM,IDISK)

#ifdef _DEBUG_
      NBUF=NPQ*JNUM
      SQFP = DNRM2_(NBUF,CHOBUF,1)
      WRITE(6,'(1X,A,I9,A,A,I2,A,A,I2,A,A,I2,A,A,F21.14)')
     &  'BATCH ',IB,   ', ',
     &  'CASE ' ,ICASE,', ',
     &  'ISYQ ' ,ISYQ, ', ',
     &  'JSYM ' ,JSYM, ', ',
     &  'DNRM2 ',SQFP
#endif
      END SUBROUTINE

************************************************************************
      SUBROUTINE CHOVEC_LOAD(CHOBUF,ICASE,ISYQ,JSYM,IB)
************************************************************************
* Read Cholesky vectors from disk.
************************************************************************
      Implicit real*8 (a-h,o-z)
#include "rasdim.fh"
#include "warnings.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "chocaspt2.fh"
      DIMENSION CHOBUF(*)

C always write the chunks to LUDRA, both for serial and parallel
      NPQ=NPQ_CHOTYPE(ICASE,ISYQ,JSYM)
      JNUM=NVLOC_CHOBATCH(IB)
      IDISK=IDLOC_CHOGROUP(ICASE,ISYQ,JSYM,IB)
      CALL DDAFILE(LUDRA,2,CHOBUF,NPQ*JNUM,IDISK)

#ifdef _DEBUG_
      NBUF=NPQ*JNUM
      SQFP = DNRM2_(NBUF,CHOBUF,1)
      WRITE(6,'(1X,A,I9,A,A,I2,A,A,I2,A,A,I2,A,A,F21.14)')
     &  'BATCH ',IB,   ', ',
     &  'CASE ' ,ICASE,', ',
     &  'ISYQ ' ,ISYQ, ', ',
     &  'JSYM ' ,JSYM, ', ',
     &  'DNRM2 ',SQFP
#endif
      END SUBROUTINE

************************************************************************
      SUBROUTINE CHOVEC_COLL(CHOBUF,ICASE,ISYQ,JSYM,IB)
************************************************************************
* Routine to gather locally available cholesky vectors and collect
* all of them on each process in case of parallel run.
************************************************************************
#ifdef _MOLCAS_MPP_
      use mpi
#endif
      IMPLICIT NONE
#include "rasdim.fh"
#include "warnings.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "chocaspt2.fh"
#include "para_info.fh"
      REAL*8 :: CHOBUF(*)
      INTEGER :: ICASE,ISYQ,JSYM,IB

#ifdef _MOLCAS_MPP_
#  include "global.fh"
#  include "mafdecls.fh"
      INTEGER*4 IERROR4,ITYPE
      INTEGER*4, PARAMETER :: ONE4 = 1
      INTEGER :: LDISP,LSIZE,LRECVBUF,LTRANSP
      INTEGER :: I,JNUM,JNUMT,NPQ,NUMSEND(1),IDISKT,IERROR
#ifdef _DEBUG_
      INTEGER :: MY_N,NOFF
      REAL*8 :: SQFP
      REAL*8, EXTERNAL :: DDOT_
#endif
#endif

#ifdef _MOLCAS_MPP_
#  ifdef _I8_
      ITYPE=MPI_INTEGER8
#  else
      ITYPE=MPI_INTEGER4
#  endif
      IF (Is_Real_Par()) THEN
C for true parallel, also communicate chunks to each process, write them
C to LUDRATOT, so first allocate memory for the fully transformed
C vectors, and for the per-process size and offset into LUDRATOT
        CALL GETMEM('DISP','ALLO','INTE',LDISP,NPROCS)
        CALL GETMEM('SIZE','ALLO','INTE',LSIZE,NPROCS)

C gather sizes of local cholesky bits
        NPQ=NPQ_CHOTYPE(ICASE,ISYQ,JSYM)
        JNUM=NVLOC_CHOBATCH(IB)
        NUMSEND(1)=NPQ*JNUM
        CALL MPI_Allgather(NUMSEND,ONE4,ITYPE,
     &       IWORK(LSIZE:LSIZE+NPROCS-1),ONE4,ITYPE,
     &       MPI_COMM_WORLD, IERROR4)
C compute offsets into the receiving array
        IWORK(LDISP)=0
        DO I=2,NPROCS
          IWORK(LDISP+I-1)=IWORK(LDISP+I-2)+IWORK(LSIZE+I-2)
        END DO

C collect the vectors
        CALL GETMEM('RECVBUF','ALLO','REAL',LRECVBUF,NFTSPC_TOT)
        CALL MPI_Barrier(MPI_COMM_WORLD, IERROR4)
        CALL MPI_Allgatherv_(CHOBUF,NUMSEND(1),MPI_REAL8,
     &       WORK(LRECVBUF),IWORK(LSIZE),IWORK(LDISP),
     &       MPI_REAL8,MPI_COMM_WORLD, IERROR)

        JNUMT=NVGLB_CHOBATCH(IB)
        ! disk offset is block offset + preceding block size
        IDISKT=IDGLB_CHOGROUP(ICASE,ISYQ,JSYM,IB)

CSVC: for RHS on demand, write transposed chovecs, else just write
        IF (RHSDIRECT) THEN
          CALL GETMEM('TRANSP','ALLO','REAL',LTRANSP,NPQ*JNUMT)
          CALL DTRANS(NPQ,JNUMT,WORK(LRECVBUF),NPQ,
     &                          WORK(LTRANSP),JNUMT)
          CALL DDAFILE(LUDRATOT,1,WORK(LTRANSP),NPQ*JNUMT,IDISKT)
          CALL GETMEM('TRANSP','FREE','REAL',LTRANSP,NPQ*JNUMT)
        ELSE
          CALL DDAFILE(LUDRATOT,1,WORK(LRECVBUF),NPQ*JNUMT,IDISKT)
        END IF

#  ifdef _DEBUG_
        WRITE(6,*) ' process block, size, offset, fingerprint'
        DO I=1,NPROCS
          MY_N = IWORK(LSIZE+I-1)
          NOFF = IWORK(LDISP+I-1)
          SQFP =DDOT_(MY_N,WORK(LRECVBUF+NOFF),1,WORK(LRECVBUF+NOFF),1)
          WRITE(6,'(A,I6,A,2I12,ES20.12)') ' [',I,'] ',MY_N,NOFF,SQFP
        END DO
#  endif
        CALL GETMEM('RECVBUF','FREE','REAL',LRECVBUF,NFTSPC_TOT)
        CALL GETMEM('DISP','FREE','INTE',LDISP,NPROCS)
        CALL GETMEM('SIZE','FREE','INTE',LSIZE,NPROCS)
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
      END SUBROUTINE

#ifdef _MOLCAS_MPP_
************************************************************************
      SUBROUTINE MPI_Allgatherv_(SENDBUF,NSEND,MPITYPES,
     &                     RCVBUF,NRCV,NOFF,MPITYPER,MPICOMM,IERROR)
************************************************************************
* Wrapper to MPI_Allgatherv dealing with ILP64 incompatibility.
************************************************************************
      use mpi
      IMPLICIT NONE
      REAL*8 SENDBUF(*), RCVBUF(*)
      INTEGER NSEND, NRCV(*),NOFF(*)

      INTEGER*4 MPITYPES, MPITYPER, MPICOMM

      INTEGER*4 NPROCS
      INTEGER*4 NSEND4
      INTEGER*4,ALLOCATABLE :: NRCV4(:),NOFF4(:)
      INTEGER*4 IERROR4
      INTEGER, PARAMETER :: I4=KIND(NSEND4)

      INTEGER :: I, IERROR
#ifdef _I8_
      INTEGER :: NRCVTOT
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
      END SUBROUTINE
#endif

      END MODULE CHOVEC_IO
