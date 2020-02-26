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
* Copyright (C) 2006, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 2006  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE DENS2_RPT2 (CI,SGM1,SGM2,G1,G2)
      IMPLICIT NONE

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "pt2_guga.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"

#include "para_info.fh"
      LOGICAL RSV_TSK

      REAL*8 CI(MXCI),SGM1(MXCI),SGM2(MXCI)
      REAL*8 G1(NLEV,NLEV),G2(NLEV,NLEV,NLEV,NLEV)

      REAL*8 GTU,GTUVX,GTUXV

      INTEGER ID
      INTEGER IST,ISU,ISV,ISX,ISTU,ISVX
      INTEGER IT,IU,IV,IX,LT,LU,LV,LX,LTU,LVX

      INTEGER ITASK,LTASK,LTASK2T,LTASK2U,NTASKS

      INTEGER ISSG,NSGM

      REAL*8, EXTERNAL :: DDOT_,DNRM2_

c Purpose: Compute the 1- and 2-electron density matrix
c arrays G1 and G2.

      CALL QENTER('DENS2_RPT2')

      CALL DCOPY_(NG1,[0.0D0],0,G1,1)
      CALL DCOPY_(NG2,[0.0D0],0,G2,1)

C For the special cases, there is no actual CI-routines involved:
c Special code for hi-spin case:
      IF(ISCF.EQ.2) THEN
        DO IT=1,NASHT
          G1(IT,IT)=1.0D00
          DO IU=1,IT-1
            G2(IT,IT,IU,IU)= 1.0D00
            G2(IU,IU,IT,IT)= 1.0D00
            G2(IT,IU,IU,IT)=-1.0D00
            G2(IU,IT,IT,IU)=-1.0D00
          END DO
        END DO
        GOTO 99
      END IF
c Special code for closed-shell:
      IF(ISCF.EQ.1 .AND. NACTEL.GT.0) THEN
        DO IT=1,NASHT
          G1(IT,IT)=2.0D00
          G2(IT,IT,IT,IT)= 2.0D00
          DO IU=1,IT-1
            G2(IT,IT,IU,IU)= 4.0D00
            G2(IU,IU,IT,IT)= 4.0D00
            G2(IT,IU,IU,IT)=-2.0D00
            G2(IU,IT,IT,IU)=-2.0D00
          END DO
        END DO
        GOTO 99
      END IF

* For the general cases, we use actual CI routine calls, and
* have to take account of orbital order.
* We will use level inices LT,LU... in these calls, but produce
* the density matrices with usual active orbital indices.
* Translation tables L2ACT and LEVEL, in pt2_guga.fh

C-SVC20100311: set up a task table with LT,LU
      nTasks=(nLev**2+nLev)/2

      CALL GETMEM ('Tasks','ALLO','INTE',lTask,2*nTasks)
      lTask2T=lTask
      lTask2U=lTask+nTasks

      iTask=0
      DO LT=1,nLev
        DO LU=1,LT
          iTask=iTask+1
          iWork(lTask2T+iTask-1)=LT
          iWork(lTask2U+iTask-1)=LU
        ENDDO
      ENDDO
      IF (iTask.NE.nTasks) WRITE(6,*) "ERROR nTasks"

      Call Init_Tsk(ID, nTasks)

C-SVC20100311: BEGIN SEPARATE TASK EXECUTION
 500  If (.NOT.Rsv_Tsk (ID,iTask)) GOTO 501

* Compute SGM1 = E_UT acting on CI, with T.ge.U,
* i.e., lowering operations. These are allowed in RAS.
C     LTU=0
C     DO 140 LT=1,NLEV
      LT=iWork(lTask2T+iTask-1)
        IST=ISM(LT)
        IT=L2ACT(LT)
C       DO 130 LU=1,LT
        LU=iWork(lTask2U+iTask-1)
C         LTU=LTU+1
          LTU=iTask
          ISU=ISM(LU)
          IU=L2ACT(LU)
          ISTU=MUL(IST,ISU)
          ISSG=MUL(ISTU,LSYM)
          NSGM=NCSF(ISSG)
C         IF(NSGM.EQ.0) GOTO 130
          IF(NSGM.EQ.0) GOTO 500
          CALL GETSGM2(LU,LT,LSYM,CI,SGM1)
          IF(ISTU.EQ.1) THEN
            GTU=DDOT_(NSGM,CI,1,SGM1,1)
            G1(IT,IU)=GTU
            G1(IU,IT)=GTU
          END IF
          LVX=0
          DO LV=1,LT
            ISV=ISM(LV)
            IV=L2ACT(LV)
            DO LX=1,LV
              LVX=LVX+1
C             IF(LVX.GT.LTU) GOTO 125
              IF(LVX.GT.LTU) GOTO 500
              ISX=ISM(LX)
              ISVX=MUL(ISV,ISX)
              IF(ISVX.NE.ISTU) GOTO 110
              IX=L2ACT(LX)
              IF(LX.EQ.LT) THEN
C then actually T=U=V=X.
                GTUVX=DDOT_(NSGM,SGM1,1,SGM1,1)
              ELSE
                CALL GETSGM2(LX,LV,ISSG,SGM1,SGM2)
                GTUVX=DDOT_(NCONF,CI,1,SGM2,1)
              END IF

              IF(LV.EQ.LX) THEN
                GTUXV=GTUVX
              ELSE
                IF(LVX.EQ.LTU) THEN
                  GTUXV=DDOT_(NSGM,SGM1,1,SGM1,1)
                ELSE
                  CALL GETSGM2(LX,LV,LSYM,CI,SGM2)
                  GTUXV=DDOT_(NSGM,SGM1,1,SGM2,1)
                END IF
              END IF
              G2(IT,IU,IV,IX)=GTUVX
              G2(IT,IU,IX,IV)=GTUXV
 110        CONTINUE
            END DO
          END DO

CSVC: The master node now continues to only handle task scheduling,
C     needed to achieve better load balancing. So it exits from the task
C     list.  It has to do it here since each process gets at least one
C     task.
#if defined (_MOLCAS_MPP_) && !defined (_GA_)
      IF (IS_REAL_PAR().AND.KING().AND.(NPROCS.GT.1)) GOTO 501
#endif

      GOTO 500
 501  CONTINUE
      CALL Free_Tsk(ID)

      CALL GETMEM ('Tasks','FREE','INTE',lTask,2*nTasks)

      CALL GAdSUM (G1,NG1)
      CALL GAdSUM (G2,NG2)

C-SVC20100311: serial part: add corrections to G2
      DO LT=1,NLEV
       IT=L2ACT(LT)
       DO LX=1,LT
        IX=L2ACT(LX)
        DO LU=LX,LT
         IU=L2ACT(LU)
         G2(IT,IU,IU,IX)=G2(IT,IU,IU,IX)-G1(IT,IX)
        END DO
       END DO
      END DO
      DO LT=2,NLEV
       IT=L2ACT(LT)
       DO LX=2,LT
        IX=L2ACT(LX)
        DO LU=1,LX-1
         IU=L2ACT(LU)
         G2(IT,IU,IU,IX)=G2(IT,IU,IU,IX)-G1(IT,IX)
        END DO
       END DO
      END DO
      LTU=0
      DO LT=1,NLEV
       IT=L2ACT(LT)
       DO LU=1,LT
        LTU=LTU+1
        IU=L2ACT(LU)
        LVX=0
        DO LV=1,LT
         IV=L2ACT(LV)
         DO LX=1,LV
          LVX=LVX+1
          IX=L2ACT(LX)
          IF(LVX.GT.LTU) GOTO 225
          GTUVX=G2(IT,IU,IV,IX)
          G2(IU,IT,IX,IV)=GTUVX
          G2(IV,IX,IT,IU)=GTUVX
          G2(IX,IV,IU,IT)=GTUVX
          GTUXV=G2(IT,IU,IX,IV)
          G2(IU,IT,IV,IX)=GTUXV
          G2(IX,IV,IT,IU)=GTUXV
          G2(IV,IX,IU,IT)=GTUXV
         END DO
        END DO
 225    CONTINUE
       END DO
      END DO
  99  CONTINUE

      IF(iPrGlb.GE.DEBUG) THEN
        WRITE(6,'("DEBUG> ",A)')
     &   "DENS2_RPT2: norms of the density matrices:"
        WRITE(6,'("DEBUG> ",A,1X,ES21.14)') "G1:", DNRM2_(NG1,G1,1)
        WRITE(6,'("DEBUG> ",A,1X,ES21.14)') "G2:", DNRM2_(NG2,G2,1)
      ENDIF

      CALL QEXIT('DENS2_RPT2')

      RETURN
      END
