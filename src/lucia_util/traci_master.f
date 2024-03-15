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
      use stdalloc, only: mma_allocate, mma_deallocate
      use  GLBBAS
      use Local_Arrays, only: Deallocate_Local_Arrays
      use rasscf_lucia, only: kvec3_length
      IMPLICIT REAL*8 (A-H,O-Z)
#include "mxpdim.fh"
#include "clunit.fh"
#include "crun.fh"
#include "cicisp.fh"
#include "cstate.fh"
#include "orbinp.fh"
#include "cands.fh"
#include "spinfo_lucia.fh"
#include "lucinp.fh"
#include "io_util.fh"
*
      DIMENSION LREC(MXNTTS),CMOMO(*)
      DIMENSION I_DUMMY(1)
      Real*8, allocatable:: VEC1(:), VEC2(:), VEC4(:)
      Real*8, Allocatable:: LCMOMO(:), LH1SAVE(:)
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
      Call mma_allocate(VEC1,LBLOCK,Label='VEC1')
      Call mma_allocate(VEC2,LBLOCK,Label='VEC2')
      Call mma_allocate(VEC3,KVEC3_LENGTH,Label='VEC3')
      Call mma_allocate(VEC4,NCONF,Label='VEC4')
*
* Transfer the CI-vector to LUC
*
      CALL BLKFO_MIN(ISSM,NREC,LREC)
      IDISK(LUC)=0
      JDISK = JOBDISK
      DO JROOT = 1,NROOT
         CALL DDAFILE(JOBIPH,2,VEC4,NCONF,JDISK)
         CALL CSDTVC(VEC4,VEC1,1,DTOC,SDREO,ISSM,0)
         IF (NTEST .GE. 50) THEN
            Write(6,*) 'CI-vector written to disk for root = ',JROOT
            CALL WRTMAT(VEC1,1,20,1,20)
            CALL XFLUSH(6)
            Write(6,*) 'Writing this to disk:'
            IOFF  = 1
            DO IREC = 1, NREC
               IF (LREC(IREC) .GE. 0) THEN
                  CALL WRTMAT(VEC1(IOFF),1,LREC(IREC),
     &                  1,LREC(IREC))
                  IOFF = IOFF + LREC(IREC)
               END IF
            END DO
            CALL XFLUSH(6)
         END IF
         CALL TODSCN(VEC1,NREC,LREC,LBLK,LUC)
         CALL ITODS([-1],1,LBLK,LUC)
      END DO

*. MO-MO transformation matrix :
      Call mma_allocate(LCMOMO,NDIM,Label='LCMOMO')
*. Copy of one-electron integrals
      Call mma_allocate(LH1SAVE,NDIM,Label='LH1SAVE')
*. We are going to mess with the one-electron integrals, take a copy
      CALL COPVEC(INT1,LH1SAVE,NDIM)
*. Set up block structure of CI space
      IATP = 1
      IBTP = 2
      CALL  Z_BLKFO(ISSPC,ISSM,IATP,IBTP,NBATCH,NBLOCK)

      Call Deallocate_Local_Arrays()
*
* The input transformation matrix contains a lot of zeros which
* is expected not to be there in Traci_Lucia, so remove them.
*
      LCMOMO(:)=0.0D0
      IOFF = 1
      IADR = 1
      ICOL = 1
      DO ISM = 1,NSMOB
         IF (NTOOBS(ISM) .GT. 0) THEN
            IROW = ICOL
            DO I = 1,NTOOBS(ISM)
               IADR = (ICOL-1)*NTOOB+IROW
               DO J = 1,NTOOBS(ISM)
                   LCMOMO(IOFF+NTOOBS(ISM)*(J-1)+I-1) =
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
        CALL COPVCD(LUC,LUSC1,VEC1,0,LBLK)
        CALL COPVCD(LUSC1,LUSC2,VEC1,1,LBLK)
*
*. Transform CI vector : Input on LUHC, output on LUDIA (!)
        CALL COPVCD(LUSC1,LUHC,VEC1,1,LBLK)
*
        CALL TRACI_LUCIA(LCMOMO,LUHC,LUDIA,ISSPC,ISSM,VEC1,VEC2)
      END DO
*     ^ End of loop over roots
      IDISK(LUDIA)=0
*
* Copy CI-vector back to MOLCAS JOBIPH file
*
      DO JROOT = 1,NROOT
         CALL FRMDSCN(VEC1,NREC,LBLK,LUDIA)
         IF (NTEST .GE. 50) THEN
            NUM_ELE = 0
            DO IREC = 1,NREC
               NUM_ELE = NUM_ELE + LREC(IREC)
            END DO
            WRITE(6,*) 'CI-Vector read from disk for root = ',JROOT
            CALL WRTMAT(VEC1,1,NUM_ELE,1,NUM_ELE)
         END IF
         CALL CSDTVC(VEC2,VEC1,2,DTOC,SDREO,ISSM,0)
         CALL DDAFILE(JOBIPH,1,VEC2,NCONF,JOBDISK)
         CALL IFRMDS(I_DUMMY,1,LBLK,LUDIA)
      END DO
      IDISK(LUDIA)=0
*
      IF(NTEST.GE.100) THEN
        DO JROOT = 1, NROOT
          CALL WRTVCD(VEC1,LUDIA,0,LBLK)
          CALL XFLUSH(6)
        END DO
      END IF
*
*. clean up time : copy 1-e integrals back in place
      CALL COPVEC(LH1SAVE,INT1,NDIM)
*
      Call mma_deallocate(VEC1)
      Call mma_deallocate(VEC2)
      Call mma_deallocate(VEC3)
      Call mma_deallocate(VEC4)
      Call mma_deallocate(LCMOMO)
      Call mma_deallocate(LH1SAVE)
*
      END
