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
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC ! Jesper

*. Output : Should outside be dimensioned as MXNTTS
      INTEGER LEN_BLK(*)
      INTEGER IDUMMY(1)
*
      IDUMMY = 0 ! jwk-cleanup
      IATP = 1
      IBTP = 2
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*. Pointers to local arrays
      IDUM=0
      NTTS = MXNTTS
      CALL GETMEM('CLBT  ','ALLO','INTE',KPCLBT ,MXNTTS)
      CALL GETMEM('CLEBT ','ALLO','INTE',KPCLEBT,MXNTTS)
      CALL GETMEM('CI1BT ','ALLO','INTE',KPCI1BT,MXNTTS)
      CALL GETMEM('CIBT  ','ALLO','INTE',KPCIBT ,8*MXNTTS)
      CALL GETMEM('CBLTP ','ALLO','INTE',KPCBLTP,NSMST)
*. Info needed for generation of block info
      CALL GETMEM('CIOIO ','ALLO','INTE',KLCIOIO,NOCTPA*NOCTPB)
      CALL IAIBCM(ISSPC,iWORK(KLCIOIO)) ! Jesper
      CALL ZBLTP(ISMOST(1,ISM),NSMST,IDC,IWORK(KPCBLTP),IDUMMY)
*. Allowed length of each batch( not important for final output )
      LBLOCK = MAX(MXSOOB,LCSBLK)
*. Batches  of C vector
      CALL PART_CIV2(IDC,IWORK(KPCBLTP),IWORK(KNSTSO(IATP)),
     &              IWORK(KNSTSO(IBTP)),
     &              NOCTPA,NOCTPB,NSMST,LBLOCK,IWORK(KLCIOIO),
     &              ISMOST(1,ISM),
     &              NBATCH,IWORK(KPCLBT),IWORK(KPCLEBT),
     &              IWORK(KPCI1BT),IWORK(KPCIBT),0,ISIMSYM)
*. Number of BLOCKS
      NBLK = IFRMR(iWORK(KPCI1BT),1,NBATCH)
     &     + IFRMR(iWORK(KPCLBT),1,NBATCH) - 1
*. Length of each block
      CALL EXTRROW(iWORK(KPCIBT),8,8,NBLK,LEN_BLK)
*
      CALL GETMEM('CLBT  ','FREE','INTE',KPCLBT ,MXNTTS)
      CALL GETMEM('CLEBT ','FREE','INTE',KPCLEBT,MXNTTS)
      CALL GETMEM('CI1BT ','FREE','INTE',KPCI1BT,MXNTTS)
      CALL GETMEM('CIBT  ','FREE','INTE',KPCIBT ,8*MXNTTS)
      CALL GETMEM('CBLTP ','FREE','INTE',KPCBLTP,NSMST)
      CALL GETMEM('CIOIO ','FREE','INTE',KLCIOIO,NOCTPA*NOCTPB)
      RETURN
      END
*
