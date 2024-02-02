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
* Copyright (C) 2001, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE BLKFO_MIN(ISM,NBLK,LEN_BLK)
      use stdalloc, only: mma_allocate, mma_deallocate
*
* Number of blocks and length of each block for CI expansion
*
* Jeppe Olsen, June 2001
*
*   Input
* =========
* ISM : Symmetry of CI expansion
*
* Output
* ======
*
* NBLK : Number of blocks in expansion
* LEN_BLK(IBLK) : Length of block IBLK
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.fh"
#include "WrkSpc.fh"
#include "cicisp.fh"
#include "stinf.fh"
#include "cstate.fh"
#include "csm.fh"
#include "strbas.fh"
#include "crun.fh"
#include "cands.fh"
      Integer, Allocatable:: CIOIO(:)
      Integer, Allocatable:: CLBT(:), CLEBT(:), CI1BT(:), CIBT(:),
     &                       CBLTP(:)

*. Output : Should outside be dimensioned as MXNTTS
      INTEGER LEN_BLK(*)
      INTEGER I_DUMMY(1)
*
      I_DUMMY(1) = 0 ! jwk-cleanup
      IATP = 1
      IBTP = 2
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*. Pointers to local arrays
      CALL mma_allocate(CLBT ,MXNTTS,Label='CLBT')
      CALL mma_allocate(CLEBT,MXNTTS,Label='CLEBT')
      CALL mma_allocate(CI1BT,MXNTTS,Label='CI1BT')
      CALL mma_allocate(CIBT ,8*MXNTTS,Label='CIBT')
      CALL mma_allocate(CBLTP,NSMST,Label='CBLTP')
*. Info needed for generation of block info
      Call mma_allocate(CIOIO,NOCTPA*NOCTPB,Label='CIOIO')
      CALL IAIBCM(ISSPC,CIOIO) ! Jesper
      CALL ZBLTP(ISMOST(1,ISM),NSMST,IDC,CBLTP,I_DUMMY)
*. Allowed length of each batch( not important for final output )
      LBLOCK = MAX(MXSOOB,LCSBLK)
*. Batches  of C vector
      CALL PART_CIV2(IDC,CBLTP,IWORK(KNSTSO(IATP)),
     &              IWORK(KNSTSO(IBTP)),
     &              NOCTPA,NOCTPB,NSMST,LBLOCK,CIOIO,
     &              ISMOST(1,ISM),
     &              NBATCH,CLBT,CLEBT,
     &              CI1BT,CIBT,0,ISIMSYM)
*. Number of BLOCKS
      NBLK = IFRMR(CI1BT,1,NBATCH)
     &     + IFRMR(CLBT,1,NBATCH) - 1
*. Length of each block
      CALL EXTRROW(CIBT,8,8,NBLK,LEN_BLK)
*
      CALL mma_deallocate(CLBT)
      CALL mma_deallocate(CLEBT)
      CALL mma_deallocate(CI1BT)
      CALL mma_deallocate(CIBT)
      CALL mma_deallocate(CBLTP)
      Call mma_deallocate(CIOIO)
      RETURN
      END
*
