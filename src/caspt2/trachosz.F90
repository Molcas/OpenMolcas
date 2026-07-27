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
! Copyright (C) Per Ake Malmqvist                                      *
!***********************************************************************

subroutine TRACHOSZ()

use Symmetry_Info, only: Mul
use CHOVEC_IO, only: IDGLB_CHOGROUP, IDLOC_CHOGROUP, NPQ_CHOTYPE, NVGLB_CHOBATCH, NVLOC_CHOBATCH, NVTOT_CHOSYM
use Para_Info, only: nProcs
use Cholesky, only: InfVec
use ChoCASPT2, only: MxCharR, MxNVC, nChSpc, nFtSpc, nHtSpc, nKsh, nPsh, NumCho_pt2
#ifdef _MOLCAS_MPP_
use ChoCASPT2, only: NFTSPC_TOT
#endif
use caspt2_global, only: do_grad, LUDRA, LUDRATOT
use general_data, only: nAsh
use caspt2_module, only: nBas, nBasT, nBtch, nBtches, nFro, nIsh, nSym
use stdalloc, only: mma_allocate, mma_MaxDBLE
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: IB, IBATCH_TOT, IBEND, IBSTA, ICASE, IDISK, ISYMA, ISYMB, ISYQ, JRED, JRED1, JRED2, JSTART, JSYM, MXFTARR, &
                     MXHTARR, MXSPC, NBATCH, NBATCH_TOT, NJSCT, NPB, NPQ, NV, NVACC, NVACT, NVECS_RED
real(kind=wp) :: Dummy(1)
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: NJSCT_TOT
#endif

! ======================================================================
! Determine sectioning size to use for the full-transformed MO vectors
! using Francesco's method.
! Cholesky vectors, and half transformed vectors, need space for
! all symmetry blocks with a specified combined symmetry.
! Fully transformed symmetry blocks are handled individually.
MXHTARR = 0
MXFTARR = 0
do JSYM=1,NSYM
  NPB = 0
  do ISYMA=1,NSYM
    ISYMB = Mul(ISYMA,JSYM)
    NPB = NPB+max(NFRO(iSymA),NISH(iSymA),NASH(iSymA))*NBAS(ISYMB)
    MXFTARR = max(MXFTARR,NPSH(ISYMA)*NKSH(ISYMB))
  end do
  MXHTARR = max(MXHTARR,NPB)
end do
MXCHARR = NBAST**2
if (do_grad) MXHTARR = MXCHARR
! MXFTARR,MXHTARR: Largest single full-transformed, half-transformed vector.
! MXCHARR: Largest possible Cholesky vector.

! What is largest possible array that can now be allocated?
call mma_MaxDBLE(MXSPC)
! Subtract 7*MXCHARR (for vector V, etc, see below).
MXSPC = MXSPC-7*MXCHARR

! Use 80% of this:
MXSPC = int(real(MXSPC,kind=wp)*0.8_wp)
! Max number of vectors that will fit in memory:
!SVC: added space for 2x the collected chovecs
MXNVC = MXSPC/(MXCHARR+MXHTARR+MXFTARR+2*nProcs*MXFTARR)
!SVC: MPI workaround: collected chovecs should not exceed 2GB
if (MXFTARR /= 0) MXNVC = min(MXNVC,2147483647/(8*nProcs*MXFTARR))
! Max number of vectors actually used in one batch:
NJSCT = 0
IBATCH_TOT = 0

do JSYM=1,NSYM
  ! Nr of batches in earlier symmetries:
  NBTCHES(JSYM) = IBATCH_TOT
  NBTCH(JSYM) = 0
  select case (NUMCHO_PT2(JSYM))
    case (0)
      NBTCH(JSYM) = 0
    case Default
      JRED1 = InfVec(1,2,jSym)
      JRED2 = InfVec(NumCho_PT2(jSym),2,jSym)
      ! Loop over the reduced sets:
      do JRED=JRED1,JRED2
        call Cho_X_nVecRS(JRED,JSYM,JSTART,NVECS_RED)
        ! It happens that a reduced set is empty:
        if (NVECS_RED == 0) cycle
        ! Reduced set JRED contains NVECS_RED vectors
        ! Reduced set JRED must be divided up into NBATCH batches
        NBATCH = 1+(NVECS_RED-1)/MXNVC
        ! Necessary number of vectors in each batch is then:
        NV = 1+(NVECS_RED-1)/NBATCH
        NJSCT = max(NV,NJSCT)
        NBTCH(JSYM) = NBTCH(JSYM)+NBATCH
      end do
  end select
  ! take maximum number of batches for this symmetry over any
  ! process, such that all procs have the same number of batches
  call GAIGOP(NBTCH(JSYM),1,'max')
  IBATCH_TOT = IBATCH_TOT+NBTCH(JSYM)
  ! Nr of batches in this symmetry:
end do

NBATCH_TOT = IBATCH_TOT

#ifdef _MOLCAS_MPP_
!SVC: take the global sum of the individual maxima
NJSCT_TOT = NJSCT
call GAIGOP_SCAL(NJSCT_TOT,'+')
#endif

! Allocate space for the Cholesky vectors:
NCHSPC = NJSCT*MXCHARR
NHTSPC = NJSCT*MXHTARR
NFTSPC = NJSCT*MXFTARR
#ifdef _MOLCAS_MPP_
NFTSPC_TOT = NJSCT_TOT*MXFTARR
#endif

#ifdef _DEBUGPRINT_
write(u6,*) ' To be allocated for ...'
write(u6,'(A,1X,I12)') '   Chol. vectors: NCHSPC     =',NCHSPC
write(u6,'(A,1X,I12)') '   half-transf  : NHTSPC     =',NHTSPC
write(u6,'(A,1X,I12)') '   full-transf:   NFTSPC     =',NFTSPC
#ifdef _MOLCAS_MPP_
write(u6,'(A,1X,I12)') '   full-transf:   NFTSPC_TOT =',NFTSPC_TOT
#endif
write(u6,*) ' Cholesky vectors per symmetry:'
write(u6,'(1X,8I12)') (NUMCHO_PT2(JSYM),JSYM=1,NSYM)
#endif

! Set up tables with the number of cholesky vectors per batch and disk
! addresses for the beginning of each batch. These arrays are accessible
! through the CHOVEC_IO module.
call MMA_ALLOCATE(NVLOC_CHOBATCH,NBATCH_TOT,Label='NVLOC_CHOBATCH')
call MMA_ALLOCATE(IDLOC_CHOGROUP,4,8,8,NBATCH_TOT,Label='IDLOC_CHOGROUP')
NVLOC_CHOBATCH = 0
IDLOC_CHOGROUP = 0

IDISK = 0
IBATCH_TOT = 0
do JSYM=1,NSYM
  if (NUMCHO_PT2(JSYM) <= 0) cycle
  JRED1 = InfVec(1,2,jSym)
  JRED2 = InfVec(NumCho_PT2(jSym),2,jSym)

  do JRED=JRED1,JRED2
    call Cho_X_nVecRS(JRED,JSYM,JSTART,NVECS_RED)
    ! It happens that a reduced set is empty:
    if (NVECS_RED == 0) cycle

    NBATCH = 1+(NVECS_RED-1)/MXNVC
    NV = 1+(NVECS_RED-1)/NBATCH
    NVACC = 0
    do IB=1,NBATCH
      IBATCH_TOT = IBATCH_TOT+1
      ! number of vectors
      NVACT = min(NVECS_RED-NVACC,NV)
      NVLOC_CHOBATCH(IBATCH_TOT) = NVACT
      NVACC = NVACC+NVACT
      ! disk address offsets
      do ISYQ=1,NSYM
        do ICASE=1,4
          NPQ = NPQ_CHOTYPE(ICASE,ISYQ,JSYM)
          IDLOC_CHOGROUP(ICASE,ISYQ,JSYM,IBATCH_TOT) = IDISK
          call DDAFILE(LUDRA,0,DUMMY,NPQ*NVACT,IDISK)
        end do
      end do
    end do
  end do
  ! set remaining batches to 0
  NBATCH = IBATCH_TOT-NBTCHES(JSYM)
  do IB=NBATCH+1,NBTCH(JSYM)
    IBATCH_TOT = IBATCH_TOT+1
    NVLOC_CHOBATCH(IBATCH_TOT) = 0
    do ISYQ=1,NSYM
      do ICASE=1,4
        IDLOC_CHOGROUP(ICASE,ISYQ,JSYM,IBATCH_TOT) = IDISK
      end do
    end do
    IDLOC_CHOGROUP(:,1:NSYM,JSYM,IBATCH_TOT) = IDISK
  end do
end do

! SVC: added workaround to get _all_ the fully transformed cholesky
! vectors onto every process. LUDRA has a counterpart LUDRATOT with
! indexing through the size NVGLB_CHOBATCH and offset IDGLB_CHOGROUP
! available from the CHOVEC_IO module.

call MMA_ALLOCATE(NVGLB_CHOBATCH,NBATCH_TOT,Label='NVGLB_CHOBATCH')
NVGLB_CHOBATCH(:) = NVLOC_CHOBATCH(:)
#ifdef _MOLCAS_MPP_
! for parrallel, sum over processes
call GAIGOP(NVGLB_CHOBATCH,NBATCH_TOT,'+')
#endif

! sum over same-symmetry batches
NVTOT_CHOSYM(:) = 0
do JSYM=1,NSYM
  IBSTA = NBTCHES(JSYM)+1
  IBEND = NBTCHES(JSYM)+NBTCH(JSYM)
  ! total size is sum over global batch sizes
  NVTOT_CHOSYM(JSYM) = sum(NVGLB_CHOBATCH(IBSTA:IBEND))
end do

call MMA_ALLOCATE(IDGLB_CHOGROUP,4,8,8,NBATCH_TOT,Label='IDGLB_CHOGROUP')

! compute offsets into all cholesky vectors
IDISK = 0
IDGLB_CHOGROUP = 0
do JSYM=1,NSYM
  IBSTA = NBTCHES(JSYM)+1
  IBEND = NBTCHES(JSYM)+NBTCH(JSYM)
  do IB=IBSTA,IBEND
    NV = NVGLB_CHOBATCH(IB)
    do ISYQ=1,NSYM
      do ICASE=1,4
        NPQ = NPQ_CHOTYPE(ICASE,ISYQ,JSYM)
        IDGLB_CHOGROUP(ICASE,ISYQ,JSYM,IB) = IDISK
        call DDAFILE(LUDRATOT,0,DUMMY,NPQ*NV,IDISK)
      end do
    end do
  end do
end do

end subroutine TRACHOSZ
