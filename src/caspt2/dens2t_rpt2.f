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
      SUBROUTINE DENS2T_RPT2 (CI1,CI2,SGM1,SGM2,G1,G2)
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par, King
#endif
      use caspt2_output, only:iPrGlb,debug
      IMPLICIT NONE

#include "rasdim.fh"
#include "caspt2.fh"
#include "pt2_guga.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"

      LOGICAL RSV_TSK

      REAL*8 CI1(MXCI),CI2(MXCI),SGM1(MXCI),SGM2(MXCI)
      REAL*8 G1(NLEV,NLEV),G2(NLEV,NLEV,NLEV,NLEV)

      REAL*8 GTU,GTUVX !! ,GTUXV

      INTEGER ID
      INTEGER IST,ISU,ISV,ISX,ISTU,ISVX
      INTEGER IT,IU,IV,IX,LT,LU,LV,LX,LVX
      integer itu,ivx

      INTEGER ITASK,LTASK,LTASK2T,LTASK2U,NTASKS

      INTEGER ISSG,NSGM

      REAL*8, EXTERNAL :: DDOT_,DNRM2_

c Purpose: Compute the 1- and 2-electron density matrix
c arrays G1 and G2.


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
      nTasks= nLev**2

      CALL GETMEM ('Tasks','ALLO','INTE',lTask,2*nTasks)
      lTask2T=lTask
      lTask2U=lTask+nTasks

      iTask=0
      DO LT=1,nLev
        DO LU=1,nLev!LT
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
          ! LTU=iTask
          ISU=ISM(LU)
          IU=L2ACT(LU)
          ISTU=MUL(IST,ISU)
          ISSG=MUL(ISTU,STSYM)
          NSGM=NCSF(ISSG)
C         IF(NSGM.EQ.0) GOTO 130
          IF(NSGM.EQ.0) GOTO 500
C         CALL GETSGM2(LT,LU,STSYM,CI1,SGM1)
C         write(6,*) "LT,LU=",lt,lu
C         do ix = 1, nsgm
C         write(6,'(i3,f20.10)') ix,sgm1(ix)
C         end do
          CALL GETSGM2(LU,LT,STSYM,CI1,SGM1)
          IF(ISTU.EQ.1) THEN
            GTU=DDOT_(NSGM,CI2,1,SGM1,1)
            G1(IT,IU)=G1(IT,IU)+GTU
C           G1(IU,IT)=GTU
          END IF
          LVX=0
          DO LV=1,NLEV!LT
            ISV=ISM(LV)
            IV=L2ACT(LV)
            DO LX=1,NLEV!LV
              LVX=LVX+1
C             IF(LVX.GT.LTU) GOTO 500
              ISX=ISM(LX)
              ISVX=MUL(ISV,ISX)
              IF(ISVX.NE.ISTU) GOTO 110
              IX=L2ACT(LX)
C             IF(LX.EQ.LT) THEN
C then actually T=U=V=X.
C               GTUVX=DDOT_(NSGM,SGM1,1,SGM1,1)
C             ELSE
                CALL GETSGM2(LX,LV,ISSG,SGM1,SGM2)
                GTUVX=DDOT_(NCONF,CI2,1,SGM2,1)
C             END IF

C             IF(LV.EQ.LX) THEN
C               GTUXV=GTUVX
C             ELSE
C               IF(LVX.EQ.LTU) THEN
C                 GTUXV=DDOT_(NSGM,SGM1,1,SGM1,1)
C               ELSE
C                 CALL GETSGM2(LX,LV,STSYM,CI,SGM2)
C                 GTUXV=DDOT_(NSGM,SGM1,1,SGM2,1)
C               END IF
C             END IF
              G2(IT,IU,IV,IX)=G2(IT,IU,IV,IX)+GTUVX
C             G2(IT,IU,IX,IV)=GTUXV
 110        CONTINUE
            END DO
          END DO

          CALL GETSGM2(LU,LT,STSYM,CI2,SGM1)
          IF(ISTU.EQ.1) THEN
            GTU=DDOT_(NSGM,CI1,1,SGM1,1)
            G1(IT,IU)=G1(IT,IU)+GTU
          END IF
          LVX=0
          DO LV=1,NLEV
            ISV=ISM(LV)
            IV=L2ACT(LV)
            DO LX=1,NLEV
              LVX=LVX+1
C             IF(LVX.GT.LTU) GOTO 500
              ISX=ISM(LX)
              ISVX=MUL(ISV,ISX)
C             IF(ISVX.NE.ISTU) GOTO 110
              IX=L2ACT(LX)
              CALL GETSGM2(LX,LV,ISSG,SGM1,SGM2)
              GTUVX=DDOT_(NCONF,CI1,1,SGM2,1)
              G2(IT,IU,IV,IX)=G2(IT,IU,IV,IX)+GTUVX
C110        CONTINUE
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

C     write(6,*) "before"
C     call sqprt(g2,nlev**2)
      Do LT = 1, NLEV
        IT = L2ACT(LT)
        Do LU = 1, NLEV
          IU = L2ACT(LU)
          Do LV = 1, NLEV
            IV = L2ACT(LV)
            G2(IT,IV,IV,IU) = G2(IT,IV,IV,IU) - G1(IT,IU)
          End Do
        End Do
      End Do
      do it = 1, nlev
      do iu = 1, nlev
      itu = it+nasht*(iu-1)
      do iv = 1, nlev
      do ix = 1, nlev
      ivx = iv+nasht*(ix-1)
       if (ivx.gt.itu) g2(iv,ix,it,iu) = g2(it,iu,iv,ix)
      end do
      end do
      end do
      end do
C-SVC20100311: serial part: add corrections to G2
C     DO LT=1,NLEV
C      IT=L2ACT(LT)
C      DO LX=1,LT
C       IX=L2ACT(LX)
C       DO LU=LX,LT
C        IU=L2ACT(LU)
C        G2(IT,IU,IU,IX)=G2(IT,IU,IU,IX)-G1(IT,IX)
C       END DO
C      END DO
C     END DO
C     DO LT=2,NLEV
C      IT=L2ACT(LT)
C      DO LX=2,LT
C       IX=L2ACT(LX)
C       DO LU=1,LX-1
C        IU=L2ACT(LU)
C        G2(IT,IU,IU,IX)=G2(IT,IU,IU,IX)-G1(IT,IX)
C       END DO
C      END DO
C     END DO
C     LTU=0
C     DO LT=1,NLEV
C      IT=L2ACT(LT)
C      DO LU=1,LT
C       LTU=LTU+1
C       IU=L2ACT(LU)
C       LVX=0
C       DO LV=1,LT
C        IV=L2ACT(LV)
C        DO LX=1,LV
C         LVX=LVX+1
C         IX=L2ACT(LX)
C         IF(LVX.GT.LTU) GOTO 225
C         GTUVX=G2(IT,IU,IV,IX)
C         G2(IU,IT,IX,IV)=GTUVX
C         G2(IV,IX,IT,IU)=GTUVX
C         G2(IX,IV,IU,IT)=GTUVX
C         GTUXV=G2(IT,IU,IX,IV)
C         G2(IU,IT,IV,IX)=GTUXV
C         G2(IX,IV,IT,IU)=GTUXV
C         G2(IV,IX,IU,IT)=GTUXV
C        END DO
C       END DO
C225    CONTINUE
C      END DO
C     END DO
  99  CONTINUE

      IF(iPrGlb.GE.DEBUG) THEN
        WRITE(6,'("DEBUG> ",A)')
     &   "DENS2_RPT2: norms of the density matrices:"
        WRITE(6,'("DEBUG> ",A,1X,ES21.14)') "G1:", DNRM2_(NG1,G1,1)
        WRITE(6,'("DEBUG> ",A,1X,ES21.14)') "G2:", DNRM2_(NG2,G2,1)
      ENDIF


      RETURN
      END
