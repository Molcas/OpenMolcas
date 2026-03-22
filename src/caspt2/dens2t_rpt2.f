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
      SUBROUTINE DENS2T_RPT2 (NLEV,NCONF,MXCI,CI1,CI2,SGM1,SGM2,G1,G2)
      use Task_Manager, only: Free_Tsk, Init_Tsk, Rsv_Tsk
      use Symmetry_Info, only: Mul
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: DEBUG
      use gugx, only: SGS, L2ACT, CIS
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: iSCF, nActEl, nAshT, STSym
      use caspt2_module, only: nG1, nG2
      use constants, only: Zero, One, Two, Four
      use definitions, only: wp, iwp, u6

      IMPLICIT NONE

      integer(kind=iwp), intent(in):: NLEV, NCONF, MXCI
      real(kind=wp), intent(in) :: CI1(NCONF), CI2(NCONF)
      real(kind=wp), intent(out) :: SGM1(MXCI), SGM2(MXCI),
     &  G1(NLEV,NLEV), G2(NLEV,NLEV,NLEV,NLEV)

      real(kind=wp) :: GTU,GTUVX !! ,GTUXV

      integer(kind=iwp) :: ID
      integer(kind=iwp) :: IST,ISU,ISV,ISX,ISTU,ISVX
      integer(kind=iwp) :: IT,IU,IV,IX,LT,LU,LV,LX,LVX
      integer(kind=iwp) :: itu,ivx
      integer(kind=iwp) :: ITASK,NTASKS
      integer(kind=iwp) :: ISSG1,ISSG2,NSGM1,NSGM2

      real(kind=wp), EXTERNAL :: DDOT_, DNRM2_
      integer(kind=iwp), ALLOCATABLE :: Task(:,:)

c Purpose: Compute the 1- and 2-electron density matrix
c arrays G1 and G2.

      G1(:,:) = Zero
      G2(:,:,:,:) = Zero

C For the special cases, there is no actual CI-routines involved:
c Special code for hi-spin case:
      IF (ISCF == 2) THEN
        DO IT = 1, NASHT
          G1(IT,IT) = One
          DO IU = 1, IT-1
            G2(IT,IT,IU,IU) =  One
            G2(IU,IU,IT,IT) =  One
            G2(IT,IU,IU,IT) = -One
            G2(IU,IT,IT,IU) = -One
          END DO
        END DO
        return
      END IF
c Special code for closed-shell:
      IF (ISCF == 1 .AND. NACTEL > 0) THEN
        DO IT = 1, NASHT
          G1(IT,IT) = Two
          G2(IT,IT,IT,IT) = Two
          DO IU = 1, IT-1
            G2(IT,IT,IU,IU) =  Four
            G2(IU,IU,IT,IT) =  Four
            G2(IT,IU,IU,IT) = -Two
            G2(IU,IT,IT,IU) = -Two
          END DO
        END DO
        return
      END IF

* For the general cases, we use actual CI routine calls, and
* have to take account of orbital order.
* We will use level inices LT,LU... in these calls, but produce
* the density matrices with usual active orbital indices.
* Translation tables L2ACT and LEVEL, in caspt2_module.F90

C-SVC20100311: set up a task table with LT,LU
      nTasks = (nLev**2+nLev)/2
      nTasks =  nLev**2

      CALL mma_allocate (Task,nTasks,2,Label='Task')

      iTask = 0
      DO LT = 1, nLev
        DO LU = 1, nLev!LT
          iTask = iTask+1
          Task(iTask,1) = LT
          Task(iTask,2) = LU
        ENDDO
      ENDDO
      IF (iTask /= nTasks) WRITE(u6,*) "ERROR nTasks"

      Call Init_Tsk(ID, nTasks)

C-SVC20100311: BEGIN SEPARATE TASK EXECUTION
      do while (Rsv_Tsk(ID,iTask))

* Compute SGM1 = E_UT acting on CI, with T.ge.U,
* i.e., lowering operations. These are allowed in RAS.
        LT=Task(iTask,1)
        IST=SGS%ISM(LT)
        IT=L2ACT(LT)
        LU=Task(iTask,2)
          ISU=SGS%ISM(LU)
          IU=L2ACT(LU)
          ISTU=Mul(IST,ISU)
          ISSG1=Mul(ISTU,STSYM)
          NSGM1=CIS%NCSF(ISSG1)
          IF (NSGM1 == 0) cycle
          CALL GETSGM2(LU,LT,STSYM,CI1,NCONF,SGM1,NSGM1)
          if (ISTU == STSYM) then
            GTU=DDOT_(NSGM1,CI2,1,SGM1,1)
            G1(IT,IU)=G1(IT,IU)+GTU
C           G1(IU,IT)=GTU
          end if
          LVX=0
          DO LV=1,NLEV!LT
            ISV=SGS%ISM(LV)
            IV=L2ACT(LV)
            DO LX=1,NLEV!LV
              LVX=LVX+1
              ISX=SGS%ISM(LX)
              ISVX=Mul(ISV,ISX)
C             IF(ISVX /= ISTU) cycle
              IX=L2ACT(LX)
              ISSG2=Mul(ISVX,ISSG1)
              NSGM2=CIS%NCSF(ISSG2)
              IF (NSGM2 == 0) cycle
              if (ISSG2 == STSYM) then
C               IF(LX.EQ.LT) THEN
C then actually T=U=V=X.
C                 GTUVX=DDOT_(NSGM,SGM1,1,SGM1,1)
C               ELSE
                  CALL GETSGM2(LX,LV,ISSG1,SGM1,NCONF,SGM2,NSGM2)
                  GTUVX=DDOT_(NSGM2,CI2,1,SGM2,1)
C               END IF

C               IF(LV.EQ.LX) THEN
C                 GTUXV=GTUVX
C               ELSE
C                 IF(LVX.EQ.LTU) THEN
C                   GTUXV=DDOT_(NSGM,SGM1,1,SGM1,1)
C                 ELSE
C                   CALL GETSGM2(LX,LV,STSYM,CI,MXCI,SGM2,NSGM)
C                   GTUXV=DDOT_(NSGM,SGM1,1,SGM2,1)
C                 END IF
C               END IF
                G2(IT,IU,IV,IX)=G2(IT,IU,IV,IX)+GTUVX
C               G2(IT,IU,IX,IV)=GTUXV
              end if
            END DO
          END DO

          CALL GETSGM2(LU,LT,STSYM,CI2,NCONF,SGM1,NSGM1)
C         IF(ISTU == 1) THEN
          if (ISTU == STSYM) then
            GTU=DDOT_(NSGM1,CI1,1,SGM1,1)
            G1(IT,IU)=G1(IT,IU)+GTU
          END IF
          LVX=0
          DO LV=1,NLEV
            ISV=SGS%ISM(LV)
            IV=L2ACT(LV)
            DO LX=1,NLEV
              LVX=LVX+1
C             IF(LVX.GT.LTU) cycle
              ISX=SGS%ISM(LX)
              IX=L2ACT(LX)
              ISVX=Mul(ISV,ISX)
C             IF (ISVX /= ISTU) cycle
              ISSG2=Mul(ISVX,ISSG1)
              NSGM2=CIS%NCSF(ISSG2)
              IF (NSGM2 == 0) cycle
              if (ISSG2 == STSYM) then
                CALL GETSGM2(LX,LV,ISSG1,SGM1,NCONF,SGM2,NSGM2)
                GTUVX=DDOT_(NSGM2,CI1,1,SGM2,1)
                G2(IT,IU,IV,IX)=G2(IT,IU,IV,IX)+GTUVX
              end if
            END DO
          END DO

CSVC: The master node now continues to only handle task scheduling,
C     needed to achieve better load balancing. So it exits from the task
C     list.  It has to do it here since each process gets at least one
C     task.

      end do

      CALL Free_Tsk(ID)

      CALL mma_deallocate(Task)

      CALL GAdGOP (G1,NG1,'+')
      CALL GAdGOP (G2,NG2,'+')

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

      IF(iPrGlb >= DEBUG) THEN
        WRITE(u6,'("DEBUG> ",A)')
     &   "DENS2_RPT2: norms of the density matrices:"
        WRITE(u6,'("DEBUG> ",A,1X,ES21.14)') "G1:", DNRM2_(NG1,G1,1)
        WRITE(u6,'("DEBUG> ",A,1X,ES21.14)') "G2:", DNRM2_(NG2,G2,1)
      ENDIF

      RETURN
      END subroutine DENS2T_RPT2
