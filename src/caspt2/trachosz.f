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
* Copyright (C) Per Ake Malmqvist                                      *
************************************************************************
      SUBROUTINE TRACHOSZ
      USE CHOVEC_IO
      USE Para_Info, ONLY: nProcs
      use ChoSwp, only: InfVec
      IMPLICIT NONE
* ----------------------------------------------------------------
#include "rasdim.fh"
#include "warnings.fh"
#include "caspt2.fh"
#include "eqsolv.fh"
#include "chocaspt2.fh"
#include "choglob.fh"
#include "WrkSpc.fh"
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif
      INTEGER IB,IBSTA,IBEND,IBATCH_TOT,NBATCH,NV
      INTEGER ICASE,ISYMA,ISYMB,ISYQ,JSYM,NPB,NPQ
      INTEGER JRED,JRED1,JRED2,JSTART
      INTEGER IDISK
      INTEGER MXFTARR,MXHTARR
      INTEGER MXSPC
      INTEGER NVACT,NVACC,NVECS_RED
**********************************************************************
*  Author : P. A. Malmqvist
**********************************************************************

* ======================================================================
* Determine sectioning size to use for the full-transformed MO vectors
* using Francesco's method.
* Cholesky vectors, and half transformed vectors, need space for
* all symmetry blocks with a specified combined symmetry.
* Fully transformed symmetry blocks are handled individually.
      MXHTARR=0
      MXFTARR=0
      DO JSYM=1,NSYM
       NPB=0
       DO ISYMA=1,NSYM
        ISYMB=MUL(ISYMA,JSYM)
        NPB=NPB+MAX(NFRO(iSymA),NISH(iSymA),NASH(iSymA))*NBAS(ISYMB)
        MXFTARR=MAX(MXFTARR,NPSH(ISYMA)*NKSH(ISYMB))
       END DO
       MXHTARR=MAX(MXHTARR,NPB)
      END DO
      MXCHARR=NBAST**2
* MXFTARR,MXHTARR: Largest single full-transformed, half-transformed vector.
* MXCHARR: Largest possible Cholesky vector.

* What is largest possible array that can now be allocated?
      CALL GETMEM('MXSPC','MAX','REAL',IP_DUMMY,MXSPC)
* Subtract 7*MXCHARR (for vector V, etc, see below).
      MXSPC=MXSPC-7*MXCHARR

* Use 80% of this:
      MXSPC=INT(DBLE(MXSPC)*0.8D0)
* Max number of vectors that will fit in memory:
CSVC: added space for 2x the collected chovecs
      MXNVC=MXSPC/(MXCHARR+MXHTARR+MXFTARR+2*nProcs*MXFTARR)
CSVC: MPI workaround: collected chovecs should not exceed 2GB
      IF (MXFTARR.NE.0) THEN
        MXNVC=MIN(MXNVC,2147483647/(8*nProcs*MXFTARR))
      END IF
* Max number of vectors actually used in one batch:
      NJSCT=0
      IBATCH_TOT=0
      DO JSYM=1,NSYM
* Nr of batches in earlier symmetries:
        NBTCHES(JSYM)=IBATCH_TOT
        NBTCH(JSYM)=0
        IF(NUMCHO_PT2(JSYM).LE.0) CYCLE
        JRED1=InfVec(1,2,jSym)
        JRED2=InfVec(NumCho_PT2(jSym),2,jSym)
* Loop over the reduced sets:
        DO JRED=JRED1,JRED2
          CALL Cho_X_nVecRS(JRED,JSYM,JSTART,NVECS_RED)
* It happens that a reduced set is empty:
          IF(NVECS_RED.eq.0) CYCLE
* Reduced set JRED contains NVECS_RED vectors
* Reduced set JRED must be divided up into NBATCH batches
          NBATCH=1+(NVECS_RED-1)/MXNVC
* Necessary number of vectors in each batch is then:
          NV=1+(NVECS_RED-1)/NBATCH
          NJSCT=MAX(NV,NJSCT)
          NBTCH(JSYM)=NBTCH(JSYM)+NBATCH
        END DO
        ! take maximum number of batches for this symmetry over any
        ! process, such that all procs have the same number of batches
        CALL GAIGOP(NBTCH(JSYM),1,'max')
        IBATCH_TOT=IBATCH_TOT+NBTCH(JSYM)
* Nr of batches in this symmetry:
      END DO
      NBATCH_TOT=IBATCH_TOT

#ifdef _MOLCAS_MPP_
CSVC: take the global sum of the individual maxima
      NJSCT_TOT=NJSCT
      CALL GAIGOP_SCAL(NJSCT_TOT,'+')
#endif

* Allocate space for the Cholesky vectors:
      NCHSPC=NJSCT*MXCHARR
      NHTSPC=NJSCT*MXHTARR
      NFTSPC=NJSCT*MXFTARR
#ifdef _MOLCAS_MPP_
      NFTSPC_TOT=NJSCT_TOT*MXFTARR
#endif

#ifdef _DEBUGPRINT_
      WRITE(6,*)' To be allocated for ...'
      WRITE(6,'(A,1X,I12)')'   Chol. vectors: NCHSPC     =',NCHSPC
      WRITE(6,'(A,1X,I12)')'   half-transf  : NHTSPC     =',NHTSPC
      WRITE(6,'(A,1X,I12)')'   full-transf:   NFTSPC     =',NFTSPC
#ifdef _MOLCAS_MPP_
      WRITE(6,'(A,1X,I12)')'   full-transf:   NFTSPC_TOT =',NFTSPC_TOT
#endif
      WRITE(6,*)' Cholesky vectors per symmetry:'
      WRITE(6,'(1X,8I12)') (NUMCHO_PT2(JSYM),JSYM=1,NSYM)
#endif

* Set up tables with the number of cholesky vectors per batch and disk
* addresses for the beginning of each batch. These arrays are accessible
* through the CHOVEC_IO module.
      ALLOCATE(NVLOC_CHOBATCH(NBATCH_TOT))
      ALLOCATE(IDLOC_CHOGROUP(4,8,8,NBATCH_TOT))
      NVLOC_CHOBATCH=0
      IDLOC_CHOGROUP=0

      IDISK=0
      IBATCH_TOT=0
      DO JSYM=1,NSYM
        IF(NUMCHO_PT2(JSYM).LE.0) CYCLE
        JRED1=InfVec(1,2,jSym)
        JRED2=InfVec(NumCho_PT2(jSym),2,jSym)

        DO JRED=JRED1,JRED2
          CALL Cho_X_nVecRS(JRED,JSYM,JSTART,NVECS_RED)
* It happens that a reduced set is empty:
          IF(NVECS_RED.eq.0) CYCLE

          NBATCH=1+(NVECS_RED-1)/MXNVC
          NV=1+(NVECS_RED-1)/NBATCH
          NVACC=0
          DO IB=1,NBATCH
            IBATCH_TOT=IBATCH_TOT+1
            ! number of vectors
            NVACT=MIN(NVECS_RED-NVACC,NV)
            NVLOC_CHOBATCH(IBATCH_TOT)=NVACT
            NVACC=NVACC+NVACT
            ! disk address offsets
            DO ISYQ=1,NSYM
              DO ICASE=1,4
                NPQ=NPQ_CHOTYPE(ICASE,ISYQ,JSYM)
                IDLOC_CHOGROUP(ICASE,ISYQ,JSYM,IBATCH_TOT)=IDISK
                CALL DDAFILE(LUDRA,0,WORK(IP_DUMMY),NPQ*NVACT,IDISK)
              END DO
            END DO
          END DO
        END DO
        ! set remaining batches to 0
        NBATCH=IBATCH_TOT-NBTCHES(JSYM)
        DO IB=NBATCH+1,NBTCH(JSYM)
          IBATCH_TOT=IBATCH_TOT+1
          NVLOC_CHOBATCH(IBATCH_TOT)=0
          DO ISYQ=1,NSYM
            DO ICASE=1,4
              IDLOC_CHOGROUP(ICASE,ISYQ,JSYM,IBATCH_TOT)=IDISK
            END DO
          END DO
        END DO
      END DO

* SVC: added workaround to get _all_ the fully transformed cholesky
* vectors onto every process. LUDRA has a counterpart LUDRATOT with
* indexing through the size NVGLB_CHOBATCH and offset IDGLB_CHOGROUP
* available from the CHOVEC_IO module.

      ALLOCATE(NVGLB_CHOBATCH(NBATCH_TOT))
      NVGLB_CHOBATCH(:)=NVLOC_CHOBATCH(:)
#ifdef _MOLCAS_MPP_
      ! for parrallel, sum over processes
      CALL GAIGOP(NVGLB_CHOBATCH,NBATCH_TOT,'+')
#endif

      ! sum over same-symmetry batches
      NVTOT_CHOSYM=0
      DO JSYM=1,NSYM
        IBSTA=NBTCHES(JSYM)+1
        IBEND=NBTCHES(JSYM)+NBTCH(JSYM)
        DO IB=IBSTA,IBEND
          ! total size is sum over global batch sizes
          NVTOT_CHOSYM(JSYM)=NVTOT_CHOSYM(JSYM)+NVGLB_CHOBATCH(IB)
        END DO
      END DO

      ALLOCATE(IDGLB_CHOGROUP(4,8,8,NBATCH_TOT))

      ! compute offsets into all cholesky vectors
      IDISK=0
      IDGLB_CHOGROUP=0
      DO JSYM=1,NSYM
        IBSTA=NBTCHES(JSYM)+1
        IBEND=NBTCHES(JSYM)+NBTCH(JSYM)
        DO IB=IBSTA,IBEND
          NV=NVGLB_CHOBATCH(IB)
          DO ISYQ=1,NSYM
            DO ICASE=1,4
              NPQ=NPQ_CHOTYPE(ICASE,ISYQ,JSYM)
              IDGLB_CHOGROUP(ICASE,ISYQ,JSYM,IB)=IDISK
              CALL DDAFILE(LUDRATOT,0,WORK(IP_DUMMY),NPQ*NV,IDISK)
            END DO
          END DO
        END DO
      END DO

      RETURN
      END

      SUBROUTINE TRACHOSZ_FREE
      USE CHOVEC_IO
      DEALLOCATE(NVLOC_CHOBATCH)
      DEALLOCATE(IDLOC_CHOGROUP)
      DEALLOCATE(NVGLB_CHOBATCH)
      DEALLOCATE(IDGLB_CHOGROUP)
      END
