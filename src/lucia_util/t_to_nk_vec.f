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
* Copyright (C) 1998, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE T_TO_NK_VEC(      T,   KORB,    ISM,   ISPC,  LUCIN,
     &                        LUCOUT,      C)
      use stdalloc, only: mma_allocate, mma_deallocate
      use Local_Arrays, only: CIBT, CBLTP, Deallocate_Local_Arrays
      use strbas
*
* Evaluate T**(NK_operator) times vector on file LUIN
* to yield vector on file LUOUT
* (NK_operator is number operator for orbital K )
*
* Note LUCIN and LUCOUT are both rewinded before read/write
* Input
* =====
*  T : Input constant
*  KORB : Orbital in symmetry order
*
*  ISM,ISPC : Symmetry and space of state on LUIN
*  C : Scratch block
*
*
* Jeppe Olsen, Feb. 98
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.fh"
#include "strinp.fh"
#include "orbinp.fh"
#include "cicisp.fh"
#include "gasstr.fh"
#include "crun.fh"
#include "csm.fh"

*. Scratch block, must hold a batch of blocks
      DIMENSION C(*)
      Integer, Allocatable:: LASTR(:), LBSTR(:)
      Integer, Allocatable:: LKAOC(:), LKBOC(:)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' T_TO_NK_VEC speaking '
        WRITE(6,*) ' ISM, ISPC = ', ISM,ISPC
      END IF
*. Set up block and batch structure of vector
      IATP = 1
      IBTP = 2
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*
      CALL Z_BLKFO(ISPC,ISM,IATP,IBTP,NBATCH,NBLOCK)
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*
      Call mma_allocate(LASTR,MXNSTR*NAEL,Label='LASTR')
      Call mma_allocate(LBSTR,MXNSTR*NBEL,Label='LBSTR')
      Call mma_allocate(LKAOC,MXNSTR,Label='LKAOC')
      Call mma_allocate(LKBOC,MXNSTR,Label='LKBOC')
*. Orbital K in type ordering
      KKORB = IREOST(KORB)
      CALL T_TO_NK_VECS   (       T,   KKORB,       C,   LUCIN,  LUCOUT,
     &                     NSTSO(IATP)%I,
     &                     NSTSO(IBTP)%I,
     &                     NBLOCK,CIBT,NAEL,NBEL,LASTR,
     &                     LBSTR,CBLTP,
     &                     NSMST,ICISTR,NTOOB,LKAOC,LKBOC)

      Call mma_deallocate(LASTR)
      Call mma_deallocate(LBSTR)
      Call mma_deallocate(LKAOC)
      Call mma_deallocate(LKBOC)

      Call Deallocate_Local_Arrays()
*
      END
