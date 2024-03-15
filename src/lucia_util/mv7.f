************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE MV7(C,HC,LUC,LUHC)
      use stdalloc, only: mma_allocate, mma_deallocate
      use strbas
*
* Outer routine for sigma vector generation
* GAS version !!!!
*
* Written in terms of RASG3/SBLOCK, May 1997
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.fh"
*
* =====
*.Input
* =====
      DIMENSION C(*),HC(*)
*
*.Definition of c and sigma
#include "cands.fh"
*
*./ORBINP/ : NACOB used
#include "orbinp.fh"
#include "cicisp.fh"
#include "cstate.fh"
#include "strinp.fh"
#include "stinf.fh"
#include "csm.fh"
#include "crun.fh"
#include "gasstr.fh"
#include "cgas.fh"
#include "lucinp.fh"
#include "cprnt.fh"
#include "oper.fh"
#include "cmxcj.fh"
      Integer, Allocatable:: SIOIO(:)
      Integer, Allocatable:: SVST(:)
! this is a the same structure as for local_arrays but it can not be
! used since lower level routines will use it.
      Integer, Allocatable:: CBLTP(:), CLBT(:), CLEBT(:), CI1BT(:),
     &                       CIBT(:)
*
      IF(ICISTR.EQ.1) THEN
        WRITE(6,*) ' MV7 does not work for ICISTR = 1'
        WRITE(6,*) ' SWITCH to ICISTR = 2,3 or program'
*        STOP ' MV7 does not work for ICISTR = 1'
        CALL SYSABENDMSG('lucia_util/mv7','Internal error',' ')
      END IF
*
      MAXK1_MX = 0
      LSCMAX_MX = 0
      IATP = 1
      IBTP = 2
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*. Arrays giving allowed type combinations
      Call mma_allocate(SIOIO,NOCTPA*NOCTPB,Label='SIOIO')
      CALL IAIBCM(ISSPC,SIOIO)
*. Arrays for additional symmetry operation
      IF(IDC.EQ.3.OR.IDC.EQ.4) THEN
        CALL mma_allocate(SVST,NSMST,Label='SVST')
        CALL SIGVST(SVST,NSMST)
      ELSE
        CALL mma_allocate(SVST,1,Label='SVST')
      END IF
*. Arrays giving block type
      Call mma_allocate(CBLTP,NSMST,Label='CBLTP')
      CALL ZBLTP(ISMOST(1,ISSM),NSMST,IDC,CBLTP,SVST)
      CALL mma_deallocate(SVST)
*. Arrays for partitioning of sigma
      NTTS = MXNTTS
      Call mma_allocate(CLBT ,NTTS,Label='CLBT')
      Call mma_allocate(CLEBT ,NTTS,Label='CLEBT')
      Call mma_allocate(CI1BT,NTTS,Label='CI1BT')
      Call mma_allocate(CIBT ,8*NTTS,Label='CIBT')
*. Batches  of C vector
c      IF (ISIMSYM.EQ.0) THEN
        LBLOCK = MXSOOB
c      ELSE
c        LBLOCK = MXSOOB_AS
c      END IF
      LBLOCK = MAX(LBLOCK,LCSBLK)
* JESPER : Should reduce I/O
      IF (ENVIRO(1:6).EQ.'RASSCF') THEN
         LBLOCK = MAX(INT(XISPSM(IREFSM,1)),MXSOOB)
         IF(PSSIGN.NE.0.0D0) LBLOCK = INT(2*XISPSM(IREFSM,1))
      ENDIF
C     WRITE(6,*) ' ISSM and ICSM in MV7 =', ISSM,ICSM
      CALL PART_CIV2(IDC,
     &               CBLTP,
     &               NSTSO(IATP)%I,
     &               NSTSO(IBTP)%I,NOCTPA,NOCTPB, NSMST,LBLOCK,
     &               SIOIO,
*
     &               ISMOST(1,ISSM),
     &               NBATCH,
     &               CLBT,
     &               CLEBT,CI1BT,CIBT,
     &               0,ISIMSYM)
      Call mma_deallocate(SIOIO)
      Call mma_deallocate(CBLTP)

      IF(ICISTR.EQ.1) THEN
       LLUC = 0
       LLUHC = 0
      ELSE
       LLUC = LUC
       LLUHC = LUHC
      END IF
*
      CALL RASSG3(C,       HC,   NBATCH,
     &            CLBT,CLEBT,
     &            CI1BT,CIBT,
     &            LLUC,  LLUHC,I_AM_OUT,N_ELIMINATED_BATCHES)
C?    WRITE(6,*) ' LSCMAX_MX = ', LSCMAX_MX
*. Eliminate local memory
      Call mma_deallocate(CLBT)
      Call mma_deallocate(CLEBT)
      Call mma_deallocate(CI1BT)
      Call mma_deallocate(CIBT)
*
      END
