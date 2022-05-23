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
      SUBROUTINE TRACI_MASTER(JOBDISK,JOBIPH,CMOMO,lrec)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "mxpdim.fh"
#include "WrkSpc.fh"
#include "clunit.fh"
#include "crun.fh"
#include "cicisp.fh"
#include "cstate.fh"
#include "glbbas.fh"
#include "orbinp.fh"
#include "cands.fh"
#include "spinfo_lucia.fh"
#include "lucinp.fh"
#include "rasscf_lucia.fh"
#include "io_util.fh"
*
      DIMENSION LREC(MXNTTS),CMOMO(*)
      DIMENSION IDUMMY(1)
*
      NTEST = 0
      LBLK  = -1
      NDIM  = NTOOB*NTOOB
      NCONF = NCSF_PER_SYM(ISSM)
* JESPER: Should reduce I/O
      LBLOCK = MAX(INT(XISPSM(IREFSM,1)),MXSOOB)
      IF(PSSIGN.NE.0.0D0) LBLOCK = INT(2.0D0*XISPSM(IREFSM,1))
*
*. The three scratch  blocks
C          GET_3BLKS(KVEC1,KVEC2,KC2)
C_REPLACED BY CALLS BELOW      CALL GET_3BLKS(KVEC1,KVEC2,KVEC3)
      CALL GETMEM('VEC1  ','ALLO','REAL',KVEC1,LBLOCK)
      CALL GETMEM('VEC2  ','ALLO','REAL',KVEC2,LBLOCK)
      CALL GETMEM('KC2   ','ALLO','REAL',KVEC3,KVEC3_LENGTH)
      CALL GETMEM('KVEC4 ','ALLO','REAL',KVEC4,NCONF)
*
* Transfer the CI-vector to LUC
*
      CALL BLKFO_MIN(ISSM,NREC,LREC)
      IDISK(LUC)=0
      JDISK = JOBDISK
      DO JROOT = 1,NROOT
         CALL DDAFILE(JOBIPH,2,WORK(KVEC4),NCONF,JDISK)
         CALL CSDTVC(WORK(KVEC4),WORK(KVEC1),1,WORK(KDTOC_POINTER),
     &                    iWORK(KSDREO_POINTER),ISSM,0)
         IF (NTEST .GE. 50) THEN
            Write(6,*) 'CI-vector written to disk for root = ',JROOT
            CALL WRTMAT(WORK(KVEC1),1,20,1,20)
            CALL XFLUSH(6)
            Write(6,*) 'Writing this to disk:'
            IOFF  = 0
            DO IREC = 1, NREC
               IF (LREC(IREC) .GE. 0) THEN
                  CALL WRTMAT(WORK(KVEC1+IOFF),1,LREC(IREC),
     &                  1,LREC(IREC))
                  IOFF = IOFF + LREC(IREC)
               END IF
            END DO
            CALL XFLUSH(6)
         END IF
         CALL TODSCN(WORK(KVEC1),NREC,LREC,LBLK,LUC)
         CALL ITODS([-1],1,LBLK,LUC)
      END DO

*. MO-MO transformation matrix :
      CALL GETMEM('CMOMO ','ALLO','REAL',KLCMOMO,NDIM)
*. Copy of one-electron integrals
      CALL GETMEM('H1SAVE','ALLO','REAL',KLH1SAVE,NDIM)
*. We are going to mess with the one-electron integrals, take a copy
      CALL COPVEC(WORK(KINT1),WORK(KLH1SAVE),NDIM)
*. Set up block structure of CI space
      IATP = 1
      IBTP = 2
      CALL  Z_BLKFO(ISSPC,ISSM,IATP,IBTP,KLCLBT,KLCLEBT,
     &      KLCI1BT,KLCIBT,KLCBLTP,NBATCH,NBLOCK)
      CALL GETMEM('CLBT  ','FREE','INTE',KLCLBT ,MXNTTS)
      CALL GETMEM('CLEBT ','FREE','INTE',KLCLEBT,MXNTTS)
      CALL GETMEM('CI1BT ','FREE','INTE',KLCI1BT,MXNTTS)
      CALL GETMEM('CIBT  ','FREE','INTE',KLCIBT ,8*MXNTTS)
      CALL GETMEM('CBLTP ','FREE','INTE',KLCBLTP,NSMST)
*
* The input transformation matrix contains a lot of zeros which
* is expected not to be there in Traci_Lucia, so remove them.
*
      CALL DCOPY_(NDIM,[0.0D0],0,WORK(KLCMOMO),1)
      IOFF = 0
      IADR = 1
      ICOL = 1
      DO ISM = 1,NSMOB
         IF (NTOOBS(ISM) .GT. 0) THEN
            IROW = ICOL
            DO I = 1,NTOOBS(ISM)
               IADR = (ICOL-1)*NTOOB+IROW
               DO J = 1,NTOOBS(ISM)
                   WORK(KLCMOMO+IOFF+NTOOBS(ISM)*(J-1)+I-1) =
     &                                       CMOMO(IADR+J-1)
               END DO
               ICOL = ICOL + 1
            END DO
            IOFF = IOFF + NTOOBS(ISM)**2
         END IF
      END DO
*
* Now the actual work
*
      IDISK(LUC)=0
      IDISK(LUDIA)=0
      DO JROOT = 1, NROOT
        IDISK(LUSC1)=0
        CALL COPVCD(LUC,LUSC1,WORK(KVEC1),0,LBLK)
        CALL COPVCD(LUSC1,LUSC2,WORK(KVEC1),1,LBLK)
*
*. Transform CI vector : Input on LUHC, output on LUDIA (!)
        CALL COPVCD(LUSC1,LUHC,WORK(KVEC1),1,LBLK)
*
        CALL TRACI_LUCIA(WORK(KLCMOMO),LUHC,LUDIA,ISSPC,ISSM,
     &             WORK(KVEC1),WORK(KVEC2))
      END DO
*     ^ End of loop over roots
      IDISK(LUDIA)=0
*
* Copy CI-vector back to MOLCAS JOBIPH file
*
      DO JROOT = 1,NROOT
         CALL FRMDSCN(WORK(KVEC1),NREC,LBLK,LUDIA)
         IF (NTEST .GE. 50) THEN
            NUM_ELE = 0
            DO IREC = 1,NREC
               NUM_ELE = NUM_ELE + LREC(IREC)
            END DO
            WRITE(6,*) 'CI-Vector read from disk for root = ',JROOT
            CALL WRTMAT(WORK(KVEC1),1,NUM_ELE,1,NUM_ELE)
         END IF
         CALL CSDTVC(WORK(KVEC2),WORK(KVEC1),2,WORK(KDTOC_POINTER),
     &               iWORK(KSDREO_POINTER),ISSM,0)
         CALL DDAFILE(JOBIPH,1,WORK(KVEC2),NCONF,JOBDISK)
         CALL IFRMDS(IDUMMY,1,LBLK,LUDIA)
      END DO
      IDISK(LUDIA)=0
*
      IF(NTEST.GE.100) THEN
        DO JROOT = 1, NROOT
          CALL WRTVCD(WORK(KVEC1),LUDIA,0,LBLK)
          CALL XFLUSH(6)
        END DO
      END IF
*
*. clean up time : copy 1-e integrals back in place
      CALL COPVEC(WORK(KLH1SAVE),WORK(KINT1),NDIM)
*
      CALL GETMEM('VEC1  ','FREE','REAL',KVEC1,LBLOCK)
      CALL GETMEM('VEC2  ','FREE','REAL',KVEC2,LBLOCK)
      CALL GETMEM('KC2   ','FREE','REAL',KVEC3,KVEC3_LENGTH)
      CALL GETMEM('KVEC4 ','FREE','REAL',KVEC4,NCONF)
      CALL GETMEM('CMOMO ','FREE','REAL',KLCMOMO,NDIM)
      CALL GETMEM('H1SAVE','FREE','REAL',KLH1SAVE,NDIM)
*
      RETURN
      END
