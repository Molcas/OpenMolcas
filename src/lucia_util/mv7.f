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
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
*
*./ORBINP/ : NACOB used
#include "orbinp.fh"
#include "cicisp.fh"
#include "strbas.fh"
#include "cstate.fh"
#include "strinp.fh"
#include "stinf.fh"
#include "csm.fh"
#include "WrkSpc.fh"
#include "crun.fh"
#include "gasstr.fh"
#include "cgas.fh"
#include "lucinp.fh"
#include "cprnt.fh"
#include "glbbas.fh"
#include "oper.fh"
      COMMON/CMXCJ/MXCJ,MAXK1_MX,LSCMAX_MX
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
      CALL GETMEM('SIOIO ','ALLO','INTE',KSIOIO,NOCTPA*NOCTPB)
      CALL IAIBCM(ISSPC,iWORK(KSIOIO))
*. Arrays for additional symmetry operation
      IF(IDC.EQ.3.OR.IDC.EQ.4) THEN
        CALL GETMEM('SVST  ','ALLO','INTE',KSVST,NSMST)
        CALL SIGVST(IWORK(KSVST),NSMST)
      ELSE
         KSVST = 1
      END IF
*. Arrays giving block type
      CALL GETMEM('SBLTP ','ALLO','INTE',KSBLTP,NSMST)
      CALL ZBLTP(ISMOST(1,ISSM),NSMST,IDC,IWORK(KSBLTP),IWORK(KSVST))
      IF(IDC.EQ.3.OR.IDC.EQ.4) THEN
        CALL GETMEM('SVST  ','FREE','INTE',KSVST,NSMST)
      END IF
*. Arrays for partitioning of sigma
      NTTS = MXNTTS
      CALL GETMEM('CLBT  ','ALLO','INTE',KLSLBT ,NTTS  )
      CALL GETMEM('CLEBT ','ALLO','INTE',KLSLEBT ,NTTS  )
      CALL GETMEM('CI1BT ','ALLO','INTE',KLSI1BT,NTTS  )
      CALL GETMEM('CIBT  ','ALLO','INTE',KLSIBT ,8*NTTS)
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
      CALL PART_CIV2(      IDC,
     &               IWORK(KSBLTP),
     &               IWORK(KNSTSO(IATP)),
     &               IWORK(KNSTSO(IBTP)),NOCTPA,NOCTPB, NSMST,LBLOCK,
     &               IWORK(KSIOIO),
*
     &               ISMOST(1,ISSM),
     &                  NBATCH,
     &               IWORK(KLSLBT),
     &               IWORK(KLSLEBT),IWORK(KLSI1BT),IWORK(KLSIBT),
     &               0,ISIMSYM)
      CALL GETMEM('SIOIO ','FREE','INTE',KSIOIO,NOCTPA*NOCTPB)
      CALL GETMEM('SBLTP ','FREE','INTE',KSBLTP,NSMST)

      IF(ICISTR.EQ.1) THEN
       LLUC = 0
       LLUHC = 0
      ELSE
       LLUC = LUC
       LLUHC = LUHC
      END IF
*
      CALL RASSG3(   C,       HC,   NBATCH,
     &            iWORK(KLSLBT),iWORK(KLSLEBT),
     &            iWORK(KLSI1BT),iWORK(KLSIBT),
     &            LLUC,  LLUHC,I_AM_OUT,N_ELIMINATED_BATCHES)
C?    WRITE(6,*) ' LSCMAX_MX = ', LSCMAX_MX
*. Eliminate local memory
      CALL GETMEM('CLBT  ','FREE','INTE',KLSLBT ,NTTS  )
      CALL GETMEM('CLEBT ','FREE','INTE',KLSLEBT ,NTTS  )
      CALL GETMEM('CI1BT ','FREE','INTE',KLSI1BT,NTTS  )
      CALL GETMEM('CIBT  ','FREE','INTE',KLSIBT ,8*NTTS)
*
*
      RETURN
      END
