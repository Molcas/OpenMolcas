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
      use strbas
      use Local_Arrays, only: CLBT, CLEBT, CI1BT, CIBT, CBLTP,
     &                        Allocate_Local_Arrays,
     &                      deallocate_Local_Arrays
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
#include "cicisp.fh"
#include "stinf.fh"
#include "cstate.fh"
#include "csm.fh"
#include "crun.fh"
#include "cands.fh"
      Integer, Allocatable:: CIOIO(:)


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
      Call Allocate_Local_Arrays(MXNTTS,NSMST)
*. Info needed for generation of block info
      Call mma_allocate(CIOIO,NOCTPA*NOCTPB,Label='CIOIO')
      CALL IAIBCM(ISSPC,CIOIO) ! Jesper
      CALL ZBLTP(ISMOST(1,ISM),NSMST,IDC,CBLTP,I_DUMMY)
*. Allowed length of each batch( not important for final output )
      LBLOCK = MAX(MXSOOB,LCSBLK)
*. Batches  of C vector
      CALL PART_CIV2(IDC,CBLTP,NSTSO(IATP)%I,
     &                         NSTSO(IBTP)%I,
     &              NOCTPA,NOCTPB,NSMST,LBLOCK,CIOIO,
     &              ISMOST(1,ISM),
     &              NBATCH,CLBT,CLEBT,
     &              CI1BT,CIBT,0,ISIMSYM)
*. Number of BLOCKS
      NBLK = IFRMR(CI1BT,1,NBATCH) + IFRMR(CLBT,1,NBATCH) - 1
*. Length of each block
      CALL EXTRROW(CIBT,8,8,NBLK,LEN_BLK)
*
      Call Deallocate_Local_Arrays()
      Call mma_deallocate(CIOIO)
      RETURN
      END
*
