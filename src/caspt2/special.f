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
      SUBROUTINE SPECIAL(G1,G2,G3,F1,F2,F3,idxG3)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION G1(NASHT,NASHT),G2(NASHT,NASHT,NASHT,NASHT),G3(*)
      DIMENSION F1(NASHT,NASHT),F2(NASHT,NASHT,NASHT,NASHT),F3(*)
      INTEGER*1 idxG3(6,*)
      INTEGER I1
      PARAMETER (I1=KIND(idxG3))
C SPECIAL-CASE ROUTINE. DELIVERS G AND F MATRICES FOR A HIGH-SPIN
C OR CLOSED-SHELL SCF CASE.
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "pt2_guga.fh"

#include "para_info.fh"
      LOGICAL RSV_TSK

      CALL QENTER('SPECIAL')

      CALL DCOPY_(NG1,[0.0D0],0,G1,1)
      CALL DCOPY_(NG2,[0.0D0],0,G2,1)
      CALL DCOPY_(NG3,[0.0D0],0,G3,1)
      CALL DCOPY_(NG1,[0.0D0],0,F1,1)
      CALL DCOPY_(NG2,[0.0D0],0,F2,1)
      CALL DCOPY_(NG3,[0.0D0],0,F3,1)

      ESUM=0.0D0
      DO I=1,NLEV
        ESUM=ESUM+ETA(I)
      END DO
C ISCF=1 for closed-shell, =2 for hispin
      OCC=2.0D0
      IF(ISCF.EQ.2) OCC=1.0D0
      DO IT=1,NASHT
        G1(IT,IT)=OCC
        LT=LEVEL(IT)
        F1(IT,IT)=(ESUM*OCC-ETA(LT))*G1(IT,IT)
      END DO

      IF(NACTEL.EQ.1) THEN
        NG3=0
        GOTO 999
      END IF

      DO IT=1,NASHT
       LT=LEVEL(IT)
       DO IU=1,NASHT
        LU=LEVEL(IU)
        G2(IT,IT,IU,IU)=G1(IT,IT)*G1(IU,IU)
        IF(IU.EQ.IT) THEN
         G2(IT,IT,IU,IU)=G2(IT,IT,IU,IU)-G1(IT,IU)
        ELSE
         G2(IT,IU,IU,IT)=-G1(IT,IT)
        END IF
        F2(IT,IT,IU,IU)=(ESUM*OCC-ETA(LT)-ETA(LU))*G2(IT,IT,IU,IU)
        F2(IT,IU,IU,IT)=(ESUM*OCC-ETA(LT)-ETA(LU))*G2(IT,IU,IU,IT)
       END DO
      END DO

      IF(NACTEL.EQ.2) THEN
        NG3=0
        GOTO 999
      END IF

      NLEV2=NLEV**2
      NLEV4=NLEV**4

      iG3=0
      nTask=NLEV4
C SVC20100908 initialize the series of tasks
      Call Init_Tsk(ID, nTask)

 500  CONTINUE
#ifdef _MOLCAS_MPP_
      IF ((NG3-iG3).LT.NLEV2) GOTO 501
#endif
 502  IF (.NOT.Rsv_Tsk(ID,iTask)) GOTO 501

      IND1=MOD(iTask-1,NLEV2)+1
      IND2=((iTask-IND1)/(NLEV2))+1
      IF(IND2.GT.IND1) GOTO 502

      IT1=MOD(IND1-1,NASHT)+1
      IU1=(IND1-IT1)/NASHT+1
      LU1=LEVEL(IU1)
      IT2=MOD(IND2-1,NASHT)+1
      IU2=(IND2-IT2)/NASHT+1
      LU2=LEVEL(IU2)

      DO IT3=1,NLEV
       DO IU3=1,NLEV
        IND3=IT3+NASHT*(IU3-1)
        IF(IND3.GT.IND2) GOTO 198
        LU3=LEVEL(IU3)
        VAL=G1(IT1,IU1)*G1(IT2,IU2)*G1(IT3,IU3)

C Here VAL is the value <PSI1|E(IT1,IU1)E(IT2,IU2)E(IT3,IU3)|PSI2>
C Add here the necessary Kronecker deltas times 2-body matrix
C elements and lower, so we get a true normal-ordered density matrix
C element.

C <PSI1|E(T1,U1,T2,U2,T3,U3)|PSI2>
C = <PSI1|E(T1,U1)E(T2,U2)E(T3,U3)|PSI2>
C -D(T3,U2)*(G2(T1,U1,T2,U3)+D(T2,U1)*G1(T1,U3))
C -D(T2,U1)*G2(T1,U2,T3,U3)
C -D(T3,U1)*G2(T2,U2,T1,U3)

        IF(IT3.EQ.IU2) THEN
          VAL=VAL-G2(IT1,IU1,IT2,IU3)
          IF(IT2.EQ.IU1) THEN
            VAL=VAL-G1(IT1,IU3)
          END IF
        END IF
        IF(IT2.EQ.IU1) THEN
          VAL=VAL-G2(IT1,IU2,IT3,IU3)
        END IF
        IF(IT3.EQ.IU1) THEN
          VAL=VAL-G2(IT2,IU2,IT1,IU3)
        END IF

C VAL is now =<PSI1|E(IT1,IU1,IT2,IU2,IT3,IU3)|PSI2>
        iG3=iG3+1
        idxG3(1,iG3)=INT(iT1,I1)
        idxG3(2,iG3)=INT(iU1,I1)
        idxG3(3,iG3)=INT(iT2,I1)
        idxG3(4,iG3)=INT(iU2,I1)
        idxG3(5,iG3)=INT(iT3,I1)
        idxG3(6,iG3)=INT(iU3,I1)
        G3(iG3)=VAL
        F3(iG3)=(ESUM*OCC-ETA(LU1)-ETA(LU2)-ETA(LU3))*VAL

 198    CONTINUE
        END DO
      END DO

CSVC: The master node now continues to only handle task scheduling,
C     needed to achieve better load balancing. So it exits from the task
C     list.  It has to do it here since each process gets at least one
C     task.
#if defined (_MOLCAS_MPP_) && !defined (_GA_)
      IF (IS_REAL_PAR().AND.KING().AND.(NPROCS.GT.1)) GOTO 501
#endif

      GO TO 500
 501  CONTINUE

C SVC2010: no more tasks, wait here for the others.
      CALL Free_Tsk(ID)

      NG3=iG3

 999  CONTINUE

      CALL QEXIT('SPECIAL')
      RETURN
      END
