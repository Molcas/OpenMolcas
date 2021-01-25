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
#include "WrkSpc.fh"
#include "strinp.fh"
#include "orbinp.fh"
#include "cicisp.fh"
#include "strbas.fh"
#include "gasstr.fh"
#include "crun.fh"
#include "csm.fh"

*. Scratch block, must hold a batch of blocks
      DIMENSION C(*)
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
      CALL Z_BLKFO(     ISPC,      ISM,     IATP,     IBTP,   KLCLBT,
     &               KLCLEBT,  KLCI1BT,   KLCIBT,  KLCBLTP,   NBATCH,
     &                NBLOCK)
C           Z_BLKFO(ISPC,ISM,IATP,IBTP,KPCLBT,KPCLEBT,
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*
      CALL GETMEM('KLASTR','ALLO','INTE',KLASTR,MXNSTR*NAEL)
      CALL GETMEM('KLBSTR','ALLO','INTE',KLBSTR,MXNSTR*NBEL)
      CALL GETMEM('KLKAOC','ALLO','INTE',KLKAOC,MXNSTR)
      CALL GETMEM('KLKBOC','ALLO','INTE',KLKBOC,MXNSTR)
*. Orbital K in type ordering
      KKORB = IREOST(KORB)
      CALL T_TO_NK_VECS   (       T,   KKORB,       C,   LUCIN,  LUCOUT,
     &                     IWORK(KNSTSO(IATP)),
     &                     IWORK(KNSTSO(IBTP)),
     &                     NBLOCK,IWORK(KLCIBT),NAEL,NBEL,IWORK(KLASTR),
     &                     IWORK(KLBSTR),IWORK(KLCBLTP),
     &                   NSMST,ICISTR,NTOOB,IWORK(KLKAOC),IWORK(KLKBOC))

      CALL GETMEM('KLASTR','FREE','INTE',KLASTR,MXNSTR*NAEL)
      CALL GETMEM('KLBSTR','FREE','INTE',KLBSTR,MXNSTR*NBEL)
      CALL GETMEM('KLKAOC','FREE','INTE',KLKAOC,MXNSTR)
      CALL GETMEM('KLKBOC','FREE','INTE',KLKBOC,MXNSTR)

      CALL GETMEM('CLBT  ','FREE','INTE',KLCLBT ,MXNTTS)
      CALL GETMEM('CLEBT ','FREE','INTE',KLCLEBT,MXNTTS)
      CALL GETMEM('CI1BT ','FREE','INTE',KLCI1BT,MXNTTS)
      CALL GETMEM('CIBT  ','FREE','INTE',KLCIBT ,8*MXNTTS)
      CALL GETMEM('CBLTP ','FREE','INTE',KLCBLTP,NSMST)
*
      RETURN
      END
