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
! Copyright (C) Steven Vancoillie                                      *
!***********************************************************************

module CHOVEC_IO

! SVC: subroutines to read/write transformed cholesky vectors from/to
! disk. These are used in tracho (where they are written) and in rhsod
! (where they are read into a global array).
!
! Cholesky vectors are distributed in batches and transformed as
! such. They are written the same way to disk:
! -----------------------------------------------------
! |IB=1               |IB=2 |... |      |IB=NBATCH_TOT|
! -----------------------------------------------------
! |Cholesky vectors   |     |    |      |   ...       |
! |(p,q,nvec_chobatch)|     |    |      |   ...       |
! |for all itk>itq    |     |    |      |   ...       |
! |and isyk,isyq=jsym |     |    |      |   ...       |
! |for jsym=1         |     |    |      | jsym=NSYM   |
! -----------------------------------------------------
! each batch contains blocks of a certain type (case and symmetry).
! Always all vectors are present, with the elements distributed over
! batches as well as processes.
!
! Cases are enumerated as follows:
!   1 =    active,inactive
!   2 =    active,active
!   3 = secondary,active
!   4 = secondary,inactive
!
! Symmetry is determined by total symmetry JSYM and the symmetry of one
! of the orbital cases ISYQ, thus all values of Mul(ISYQ,JSYM).
!
! When reading the vectors, we group all batches of the same type
! together, so that we have all vectors of one type as one block,
! as this is how they are used to compute the integrals for RHS.

! indexing of distributed cholesky vectors
!-----------------------------------------
! The cholesky vectors are distributed over batches and processed
! like that by tracho2. A batch has a specified number of cholesky
! vectors (batches for different JSYM are different, so are the
! number of vectors in them). The offset of a batch on disk is
! then determined by the type (case,isyq,jsym) and total batch
! number (ib).
! NVLOC_CHOBATCH, IDLOC_CHOGROUP

! indexing of collected cholesky vectors
!---------------------------------------
! The cholesky vectors can be collected together on disk from
! different processes. The combined sizes are stored in separate
! arrays.
! NVGLB_CHOBATCH, IDGLB_CHOGROUP

! total amount of cholesky vectors in a certain symmetry
! NVTOT_CHOSYM

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: NVTOT_CHOSYM(8)
integer(kind=iwp), allocatable :: IDGLB_CHOGROUP(:,:,:,:), IDLOC_CHOGROUP(:,:,:,:), NVGLB_CHOBATCH(:), NVLOC_CHOBATCH(:)

public :: CHOVEC_COLL, CHOVEC_LOAD, CHOVEC_READ, CHOVEC_SAVE, CHOVEC_SIZE, IDGLB_CHOGROUP, IDLOC_CHOGROUP, NPQ_CHOTYPE, &
          NVGLB_CHOBATCH, NVLOC_CHOBATCH, NVTOT_CHOSYM

contains

function NPQ_CHOTYPE(ICASE,ISYQ,JSYM)
  !*********************************************************************
  ! Compute the number of orbital pairs for a given case (valid pair of
  ! inactive,active,secondary), total symmetry JSYM and component symmetry
  ! ISYQ. Note that the component Q which determines the symmetry block is
  ! the _lower_ orbital partition (e.g. inactive for active,inactive),
  ! which is also the slowest varying index of the pair P,Q.
  !*********************************************************************

  use Symmetry_Info, only: Mul
  use caspt2_module, only: nAsh, nIsh, nSSh

  integer(kind=iwp) :: NPQ_CHOTYPE
  integer(kind=iwp), intent(in) :: ICASE, ISYQ, JSYM
  integer(kind=iwp) :: ISYP, NP, NQ

  ISYP = Mul(ISYQ,JSYM)
  select case (ICASE)
    case (1)
      NP = NASH(ISYP)
      NQ = NISH(ISYQ)
    case (2)
      NP = NASH(ISYP)
      NQ = NASH(ISYQ)
    case (3)
      NP = NSSH(ISYP)
      NQ = NASH(ISYQ)
    case (4)
      NP = NSSH(ISYP)
      NQ = NISH(ISYQ)
    case default
      NP = 0
      NQ = 0
      call SYSABENDMSG('NPQ_CHOTYPE','invalid case number','')
  end select
  NPQ_CHOTYPE = NP*NQ

end function NPQ_CHOTYPE

subroutine CHOVEC_SIZE(ICASE,NCHOBUF,IOFF)
  !*********************************************************************
  ! Allocate a buffer to hold all cholesky vectors of type ITK,ITQ
  !*********************************************************************

  use Symmetry_Info, only: Mul
  use caspt2_module, only: nSym

  integer(kind=iwp), intent(in) :: ICASE
  integer(kind=iwp), intent(out) :: NCHOBUF, IOFF(8,8)
  integer(kind=iwp) :: ISYK, ISYQ, JSYM, NPQ, NVTOT

  NCHOBUF = 0
  do JSYM=1,NSYM
    NVTOT = NVTOT_CHOSYM(JSYM)
    do ISYQ=1,NSYM
      ISYK = Mul(ISYQ,JSYM)
      IOFF(ISYK,ISYQ) = NCHOBUF
      NPQ = NPQ_CHOTYPE(ICASE,ISYQ,JSYM)
      NCHOBUF = NCHOBUF+NPQ*NVTOT
    end do
  end do

end subroutine CHOVEC_SIZE

subroutine CHOVEC_READ(ICASE,CHOBUF,nCHOBUF)
  !*********************************************************************
  ! Read (transposed) cholesky vectors from disk, they
  ! are indexed as CHOBUF(IVEC,IQ,IK)
  !*********************************************************************

  use caspt2_module, only: nSym, nBtches, nBtch
  use caspt2_global, only: LUDRA
# ifdef _MOLCAS_MPP_
  use Para_Info, only: Is_Real_Par
  use caspt2_global, only: LUDRATOT
# endif
  use stdalloc, only: mma_allocate, mma_deallocate

  integer(kind=iwp), intent(in) :: ICASE, nCHOBUF
  real(kind=wp), intent(inout) :: CHOBUF(nCHOBUF)
  integer(kind=iwp) :: IB, IBEND, IBOFF, IBSTA, IDISK, IOFF, ISYQ, J, JSYM, NBUF, NPQ, NV, NVTOT
  real(kind=wp), allocatable :: BUF(:)
# ifdef _MOLCAS_MPP_
# include "global.fh"
# include "mafdecls.fh"
# endif

  IOFF = 0
  do JSYM=1,NSYM
    NVTOT = NVTOT_CHOSYM(JSYM)
    IBSTA = NBTCHES(JSYM)+1
    IBEND = NBTCHES(JSYM)+NBTCH(JSYM)
    do ISYQ=1,NSYM
      NPQ = NPQ_CHOTYPE(ICASE,ISYQ,JSYM)
      IBOFF = 0
      do IB=IBSTA,IBEND
        NV = NVGLB_CHOBATCH(IB)
        NBUF = NPQ*NV
        call mma_allocate(BUF,NBUF,LABEL='BUF')
        IDISK = IDGLB_CHOGROUP(ICASE,ISYQ,JSYM,IB)
#       ifdef _MOLCAS_MPP_
        if (Is_Real_Par()) then
          ! cholesky vectors already transposed
          call DDAFILE(LUDRATOT,2,BUF,NBUF,IDISK)
          do J=1,NPQ
            CHOBUF(IOFF+IBOFF+NVTOT*(J-1)+1:IOFF+IBOFF+NVTOT*(J-1)+NV) = BUF(NV*(J-1)+1:NV*J)
          end do
        else
#       endif
          ! cholesky vectors not transposed
          call DDAFILE(LUDRA,2,BUF,NBUF,IDISK)
          do J=1,NPQ
            CHOBUF(IOFF+IBOFF+NVTOT*(J-1)+1:IOFF+IBOFF+NVTOT*(J-1)+NV) = BUF(J:J+NPQ*(NV-1):NPQ)
          end do
#       ifdef _MOLCAS_MPP_
        end if
#       endif
        call mma_deallocate(BUF)
        IBOFF = IBOFF+NV
      end do
      IOFF = IOFF+NVTOT*NPQ
    end do
  end do

end subroutine CHOVEC_READ

subroutine CHOVEC_SAVE(CHOBUF,NCHOBUF,ICASE,ISYQ,JSYM,IB)
  !*********************************************************************
  ! Write Cholesky vectors to disk.
  !*********************************************************************

  use caspt2_global, only: LUDRA
# ifdef _DEBUGPRINT_
  use Definitions, only: u6
# endif

  integer(kind=iwp), intent(in) :: NCHOBUF, ICASE, ISYQ, JSYM, IB
  real(kind=wp), intent(inout) :: CHOBUF(NCHOBUF)
  integer(kind=iwp) :: IDISK, JNUM, NPQ
# ifdef _DEBUGPRINT_
  integer(kind=iwp) :: NBUF
  real(kind=wp) :: SQFP
  real(kind=wp), external :: DNRM2_
# endif

  ! always write the chunks to LUDRA, both for serial and parallel
  NPQ = NPQ_CHOTYPE(ICASE,ISYQ,JSYM)
  JNUM = NVLOC_CHOBATCH(IB)
  IDISK = IDLOC_CHOGROUP(ICASE,ISYQ,JSYM,IB)
  call DDAFILE(LUDRA,1,CHOBUF,NPQ*JNUM,IDISK)

# ifdef _DEBUGPRINT_
  NBUF = NPQ*JNUM
  SQFP = DNRM2_(NBUF,CHOBUF,1)
  write(u6,'(1X,A,I9,A,A,I2,A,A,I2,A,A,I2,A,A,F21.14)') &
    'BATCH ',IB,', ','CASE ',ICASE,', ','ISYQ ',ISYQ,', ','JSYM ',JSYM,', ','DNRM2 ',SQFP
# endif

end subroutine CHOVEC_SAVE

subroutine CHOVEC_LOAD(CHOBUF,NCHOBUF,ICASE,ISYQ,JSYM,IB)
  !*********************************************************************
  ! Read Cholesky vectors from disk.
  !*********************************************************************

  use caspt2_global, only: LUDRA
# ifdef _DEBUGPRINT_
  use Definitions, only: u6
# endif

  integer(kind=iwp), intent(in) :: NCHOBUF, ICASE, ISYQ, JSYM, IB
  real(kind=wp), intent(out) :: CHOBUF(NCHOBUF)
  integer(kind=iwp) :: IDISK, JNUM, NPQ
# ifdef _DEBUGPRINT_
  integer(kind=iwp) :: NBUF
  real(kind=wp) :: SQFP
  real(kind=wp), external :: DNRM2_
# endif

  ! always write the chunks to LUDRA, both for serial and parallel
  NPQ = NPQ_CHOTYPE(ICASE,ISYQ,JSYM)
  JNUM = NVLOC_CHOBATCH(IB)
  IDISK = IDLOC_CHOGROUP(ICASE,ISYQ,JSYM,IB)
  call DDAFILE(LUDRA,2,CHOBUF,NPQ*JNUM,IDISK)

# ifdef _DEBUGPRINT_
  NBUF = NPQ*JNUM
  SQFP = DNRM2_(NBUF,CHOBUF,1)
  write(u6,'(1X,A,I9,A,A,I2,A,A,I2,A,A,I2,A,A,F21.14)') &
    'BATCH ',IB,', ','CASE ',ICASE,', ','ISYQ ',ISYQ,', ','JSYM ',JSYM,', ','DNRM2 ',SQFP
# endif

end subroutine CHOVEC_LOAD

subroutine CHOVEC_COLL(CHOBUF,NCHOBUF,ICASE,ISYQ,JSYM,IB)
  !*********************************************************************
  ! Routine to gather locally available cholesky vectors and collect
  ! all of them on each process in case of parallel run.
  !*********************************************************************

# ifdef _MOLCAS_MPP_
  use MPI_Wrapper, only: MPI_COMM_WORLD, MPI_INTEGER, MPI_REAL8
  use Para_Info, only: Is_Real_Par, nProcs
  use caspt2_global, only: LUDRATOT
  use caspt2_module, only: RHSDirect
  use chocaspt2, only: NFTSPC_TOT
  use stdalloc, only: mma_allocate, mma_deallocate
  use Definitions, only: MPIInt
# endif

  integer(kind=iwp), intent(in) :: NCHOBUF, ICASE, ISYQ, JSYM, IB
  real(kind=wp), intent(inout) :: CHOBUF(NCHOBUF)
# ifdef _MOLCAS_MPP_
  integer(kind=iwp) :: I, IDISKT, IERROR, JNUM, JNUMT, NPQ, NUMSEND(1)
  integer(kind=MPIInt) :: IERROR4, ITYPE
  integer(kind=MPIInt), allocatable :: DISP(:), SIZ(:)
  real(kind=wp), allocatable :: RECVBUF(:), TRANSP(:)
# ifdef _DEBUGPRINT_
  integer(kind=iwp) :: MY_N, NOFF
  real(kind=wp) :: SQFP
  real(kind=wp), external :: DDOT_
# endif
# include "global.fh"
# include "mafdecls.fh"
# endif

# ifdef _MOLCAS_MPP_
  ITYPE = MPI_INTEGER
  if (Is_Real_Par()) then
    ! for true parallel, also communicate chunks to each process, write them
    ! to LUDRATOT, so first allocate memory for the fully transformed
    ! vectors, and for the per-process size and offset into LUDRATOT
    call mma_allocate(DISP,NPROCS,Label='DISP')
    call mma_allocate(SIZ,NPROCS,LABEL='SIZ')

    ! gather sizes of local cholesky bits
    NPQ = NPQ_CHOTYPE(ICASE,ISYQ,JSYM)
    JNUM = NVLOC_CHOBATCH(IB)
    NUMSEND(1) = NPQ*JNUM
    call MPI_Allgather(NUMSEND,1_MPIInt,ITYPE,SIZ(1:NPROCS),1_MPIInt,ITYPE,MPI_COMM_WORLD,IERROR4)
    ! compute offsets into the receiving array
    DISP(1) = 0
    do I=2,NPROCS
      DISP(I) = DISP(I-1)+SIZ(I-1)
    end do

    ! collect the vectors
    call mma_allocate(RECVBUF,NFTSPC_TOT,Label='RECVBUF')
    call MPI_Barrier(MPI_COMM_WORLD,IERROR4)
    call MPI_Allgatherv_(CHOBUF,NCHOBUF,NUMSEND(1),MPI_REAL8,RECVBUF,NFTSPC_TOT,SIZ,DISP,NPROCS,MPI_REAL8,MPI_COMM_WORLD,IERROR)

    JNUMT = NVGLB_CHOBATCH(IB)
    ! disk offset is block offset + preceding block size
    IDISKT = IDGLB_CHOGROUP(ICASE,ISYQ,JSYM,IB)

    !SVC: for RHS on demand, write transposed chovecs, else just write
    if (RHSDIRECT) then
      call mma_allocate(TRANSP,NPQ*JNUMT,Label='TRANSP')
      call DTRANS(NPQ,JNUMT,RECVBUF,NPQ,TRANSP,JNUMT)
      call DDAFILE(LUDRATOT,1,TRANSP,NPQ*JNUMT,IDISKT)
      call mma_deallocate(TRANSP)
    else
      call DDAFILE(LUDRATOT,1,RECVBUF,NPQ*JNUMT,IDISKT)
    end if

#   ifdef _DEBUGPRINT_
    write(u6,*) ' process block, size, offset, fingerprint'
    do I=1,NPROCS
      MY_N = SIZ(I)
      NOFF = 1+DISP(I)
      SQFP = DDOT_(MY_N,RECVBUF(NOFF:),1,RECVBUF(NOFF:),1)
      write(u6,'(A,I6,A,2I12,ES20.12)') ' [',I,'] ',MY_N,NOFF,SQFP
    end do
#   endif
    call mma_deallocate(RECVBUF)
    call mma_deallocate(DISP)
    call mma_deallocate(SIZ)
  end if
# else
# include "macros.fh"
  unused_var(CHOBUF)
  unused_var(ICASE)
  unused_var(ISYQ)
  unused_var(JSYM)
  unused_var(IB)
# endif

end subroutine CHOVEC_COLL

#ifdef _MOLCAS_MPP_
subroutine MPI_Allgatherv_(SENDBUF,NSENDBUF,NSEND,MPITYPES,RCVBUF,NRCVBUF,NRCV,NOFF,MPROCS,MPITYPER,MPICOMM,IERROR)
  !*********************************************************************
  ! Wrapper to MPI_Allgatherv dealing with ILP64 incompatibility.
  !*********************************************************************

  use MPI_Wrapper, only: MPI_COMM_WORLD
  use stdalloc, only: mma_allocate, mma_deallocate
  use Definitions, only: MPIInt
# ifdef _I8_
  use Definitions, only: u6
# endif

  integer(kind=iwp), intent(in) :: NSENDBUF, NSEND, NRCVBUF, MPROCS
  real(kind=wp), intent(inout) :: SENDBUF(NSENDBUF), RCVBUF(NRCVBUF)
  integer(kind=MPIInt), intent(in) :: MPITYPES, NRCV(MPROCS), NOFF(MPROCS), MPITYPER, MPICOMM
  integer(kind=iwp), intent(out) :: IERROR
  integer(kind=MPIInt) :: IERROR4, NPROCS, NSEND4
  integer(kind=MPIInt), allocatable :: NRCV4(:), NOFF4(:)
# ifdef _I8_
  integer(kind=iwp) :: NRCVTOT
# endif

  call MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS,IERROR4)

# ifdef _I8_
  NRCVTOT = sum(NRCV(1:NPROCS))
  if (8*NRCVTOT > 2147483647) then
    write(u6,'(1X,A)') 'MPI_Allgatherv: total rcv buf > 2**31-1'
    write(u6,'(1X,A)') 'workaround to avoid buffers >2GB failed'
    write(u6,'(1X,A)') '-> please report this as a bug!'
    call ABEND()
  end if
# endif

  call MMA_ALLOCATE(NRCV4,int(NPROCS,kind=iwp),Label='NRCV4')
  call MMA_ALLOCATE(NOFF4,int(NPROCS,kind=iwp),Label='NOFF4')
  NSEND4 = int(NSEND,kind=MPIInt)
  NRCV4(:) = int(NRCV(1:NPROCS),kind=MPIInt)
  NOFF4(:) = int(NOFF(1:NPROCS),kind=MPIInt)
  call MPI_Allgatherv(SENDBUF,NSEND4,MPITYPES,RCVBUF,NRCV4,NOFF4,MPITYPER,MPICOMM,IERROR4)

  IERROR = IERROR4
  call MMA_DEALLOCATE(NRCV4)
  call MMA_DEALLOCATE(NOFF4)

end subroutine MPI_Allgatherv_
#endif

end module CHOVEC_IO
