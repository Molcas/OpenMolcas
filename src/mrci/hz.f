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
      SUBROUTINE HZ(ARR)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "SysDef.fh"

#include "mrci.fh"
      PARAMETER (IX1F=1,IX2F=2,IRR=3,IX1R=4,IX2R=5,IX1X1=6,IX2X1=7,
     *           IX2X2=8,IFDF=9,IFDR=10,IRDR=11)
      DIMENSION TMP(MXVEC,MXVEC)
      DIMENSION ARR(NRROOT,NRROOT,11)
*
* THIS SUBROUTINE FORMS THE OVERLAP AND H-ZERO MATRIX ELEMENTS
* IN THE BASIS OF PSI, RHO, XI1, AND XI2 FUNCTIONS.
*
C     WRITE(6,*)
C     WRITE(6,*)' CHECK PRINTS IN HZ.'
C     WRITE(6,*)' X1F ARRAY:'
C     DO 1001 I=1,NRROOT
C       WRITE(6,'(1X,5F15.6)')(ARR(I,J,IX1F),J=1,NRROOT)
C1001 CONTINUE
C     WRITE(6,*)' X2F ARRAY:'
C     DO 1002 I=1,NRROOT
C       WRITE(6,'(1X,5F15.6)')(ARR(I,J,IX2F),J=1,NRROOT)
C1002 CONTINUE
C     WRITE(6,*)' RR ARRAY:'
C     DO 1003 I=1,NRROOT
C       WRITE(6,'(1X,5F15.6)')(ARR(I,J,IRR ),J=1,NRROOT)
C1003 CONTINUE
C     WRITE(6,*)' X1R ARRAY:'
C     DO 1004 I=1,NRROOT
C       WRITE(6,'(1X,5F15.6)')(ARR(I,J,IX1R),J=1,NRROOT)
C1004 CONTINUE
C     WRITE(6,*)' X2R ARRAY:'
C     DO 1005 I=1,NRROOT
C       WRITE(6,'(1X,5F15.6)')(ARR(I,J,IX2R),J=1,NRROOT)
C1005 CONTINUE
C     WRITE(6,*)' X1X1 ARRAY:'
C     DO 1006 I=1,NRROOT
C       WRITE(6,'(1X,5F15.6)')(ARR(I,J,IX1X1),J=1,NRROOT)
C1006 CONTINUE
C     WRITE(6,*)' X2X1 ARRAY:'
C     DO 1007 I=1,NRROOT
C       WRITE(6,'(1X,5F15.6)')(ARR(I,J,IX2X1),J=1,NRROOT)
C1007 CONTINUE
C     WRITE(6,*)' X2X2 ARRAY:'
C     DO 1008 I=1,NRROOT
C       WRITE(6,'(1X,5F15.6)')(ARR(I,J,IX2X2),J=1,NRROOT)
C1008 CONTINUE
C     WRITE(6,*)' FDF ARRAY:'
C     DO 1009 I=1,NRROOT
C       WRITE(6,'(1X,5F15.6)')(ARR(I,J,IFDF),J=1,NRROOT)
C1009 CONTINUE
C     WRITE(6,*)' FDR ARRAY:'
C     DO 1010 I=1,NRROOT
C       WRITE(6,'(1X,5F15.6)')(ARR(I,J,IFDR),J=1,NRROOT)
C1010 CONTINUE
C     WRITE(6,*)' RDR ARRAY:'
C     DO 1011 I=1,NRROOT
C       WRITE(6,'(1X,5F15.6)')(ARR(I,J,IRDR),J=1,NRROOT)
C1011 CONTINUE
C FIRST, CREATE OVERLAP MATRIX, AND INITIALIZE HZERO MATRIX WITH ALL
C TERMS THAT DO NOT REQUIRE MATRIX MULTIPLIES:
      DO 10 I1=1,NRROOT
        I2=I1+NRROOT
        I3=I1+2*NRROOT
        I4=I1+3*NRROOT
        DO 11 J1=1,NRROOT
          J2=J1+NRROOT
          J3=J1+2*NRROOT
          J4=J1+3*NRROOT
          SZERO(I1,J1)=0.0D00
          SZERO(I2,J1)=0.0D00
          SZERO(I3,J1)=ARR(I1,J1,IX1F)
          SZERO(I4,J1)=ARR(I1,J1,IX2F)
          SZERO(I2,J2)=ARR(I1,J1,IRR)
          SZERO(I3,J2)=ARR(I1,J1,IX1R)
          SZERO(I4,J2)=ARR(I1,J1,IX2R)
          SZERO(I3,J3)=ARR(I1,J1,IX1X1)
          SZERO(I4,J3)=ARR(I1,J1,IX2X1)
          SZERO(I4,J4)=ARR(I1,J1,IX2X2)
          HZERO(I1,J1)=0.0D00
          IF(I1.EQ.J1) THEN
            SZERO(I1,J1)=1.0D00
            HZERO(I1,J1)=ESMALL(I1)
          END IF
          HZERO(I2,J1)=ARR(I1,J1,IRR)
          HZERO(I3,J1)=ARR(I1,J1,IX1F)*ESMALL(J1)+ARR(I1,J1,IX1R)
          HZERO(I4,J1)=ARR(I1,J1,IX2F)*ESMALL(J1)+ARR(I1,J1,IX2R)
          HZERO(I2,J2)=ARR(I1,J1,IRDR)
          HZERO(I3,J2)=ESMALL(I1)*ARR(I1,J1,IX1R)
          HZERO(I4,J2)=ESMALL(I1)*ARR(I1,J1,IX2R)+ARR(I1,J1,IRR)
          HZERO(I3,J3)=ARR(I1,J1,IX1X1)*ESMALL(J1)-ARR(J1,I1,IX1F)
          HZERO(I4,J3)=ARR(I1,J1,IX2X1)*ESMALL(J1)
          HZERO(I4,J4)=ARR(I1,J1,IX2X2)*ESMALL(J1)+ARR(I1,J1,IX2R)
11      CONTINUE
10    CONTINUE
      IO1=0
      IO2=NRROOT
      IO3=2*NRROOT
      IO4=3*NRROOT
      DO 20 I=1,NRROOT
        DO 22 J=1,NRROOT
          SUM1=HZERO(IO3+I,IO2+J)
          SUM2=HZERO(IO4+I,IO2+J)
          DO 21 K=1,NRROOT
            SUM1=SUM1+ARR(I,K,IX1F)*(ARR(K,J,IRR)-ARR(K,J,IFDR))
            SUM2=SUM2+ARR(I,K,IX2F)*(ARR(K,J,IRR)-ARR(K,J,IFDR))
21        CONTINUE
          HZERO(IO3+I,IO2+J)=SUM1
          HZERO(IO4+I,IO2+J)=SUM2
22      CONTINUE
20    CONTINUE
      DO 30 I=1,NRROOT
        DO 32 J=1,NRROOT
          SUM1=0.0D00
          SUM2=0.0D00
          DO 31 K=1,NRROOT
            SUM1=SUM1+ARR(I,K,IX1F)*ARR(J,K,IX1F)
            SUM2=SUM2+ARR(I,K,IX2F)*ARR(J,K,IX1F)
31        CONTINUE
          HZERO(IO3+I,IO3+J)=HZERO(IO3+I,IO3+J)-SUM1*ESMALL(J)
          HZERO(IO3+J,IO3+I)=HZERO(IO3+J,IO3+I)-SUM1*ESMALL(J)
          HZERO(IO4+I,IO3+J)=HZERO(IO4+I,IO3+J)-SUM2*ESMALL(J)
32      CONTINUE
30    CONTINUE
      DO 40 I=1,NRROOT
        DO 42 J=1,NRROOT
          SUM1=0.0D00
          SUM2=0.0D00
          DO 41 K=1,NRROOT
            SUM1=SUM1+ARR(I,K,IX2F)*ARR(J,K,IX1F)
            SUM2=SUM2+ARR(I,K,IX2F)*ARR(J,K,IX2F)
41        CONTINUE
          HZERO(IO4+I,IO3+J)=HZERO(IO4+I,IO3+J)-SUM1*ESMALL(I)
          HZERO(IO4+I,IO4+J)=HZERO(IO4+I,IO4+J)-SUM2*ESMALL(I)
          HZERO(IO4+J,IO4+I)=HZERO(IO4+J,IO4+I)-SUM2*ESMALL(I)
42      CONTINUE
40    CONTINUE
      DO 50 I=1,NRROOT
        DO 52 J=1,NRROOT
          SUM=0.0D00
          DO 51 K=1,NRROOT
            SUM=SUM+ARR(I,K,IFDF)*ARR(J,K,IX1F)
51        CONTINUE
          TMP(I,J)=SUM
52      CONTINUE
50    CONTINUE
      DO 60 I=1,NRROOT
        DO 62 J=1,NRROOT
          SUM1=HZERO(IO3+I,IO3+J)
          SUM2=HZERO(IO4+I,IO3+J)
          DO 61 K=1,NRROOT
            SUM1=SUM1+ARR(I,K,IX1F)*TMP(K,J)
            SUM2=SUM2+ARR(I,K,IX2F)*TMP(K,J)
61        CONTINUE
          HZERO(IO3+I,IO3+J)=SUM1
          HZERO(IO4+I,IO3+J)=SUM2
62      CONTINUE
60    CONTINUE
      DO 70 I=1,NRROOT
        DO 72 J=1,NRROOT
          SUM=0.0D00
          DO 71 K=1,NRROOT
            SUM=SUM+ARR(I,K,IFDF)*ARR(J,K,IX2F)
71        CONTINUE
          TMP(I,J)=SUM
72      CONTINUE
70    CONTINUE
      DO 80 I=1,NRROOT
        DO 82 J=1,NRROOT
          SUM2=HZERO(IO4+I,IO4+J)
          DO 81 K=1,NRROOT
            SUM2=SUM2+ARR(I,K,IX2F)*TMP(K,J)
81        CONTINUE
          HZERO(IO4+I,IO4+J)=SUM2
82      CONTINUE
80    CONTINUE
      DO 90 I=1,NRROOT
        DO 92 J=1,NRROOT
          TMP(I,J)=ESMALL(I)*ARR(J,I,IX1F)+ARR(J,I,IX1R)
92      CONTINUE
90    CONTINUE
      DO 100 I=1,NRROOT
        DO 102 J=1,NRROOT
          SUM1=HZERO(IO3+I,IO3+J)
          SUM2=HZERO(IO4+I,IO3+J)
          DO 101 K=1,NRROOT
            SUM1=SUM1+ARR(I,K,IX1F)*TMP(K,J)+ARR(I,K,IX1R)*ARR(J,K,IX1F)
            SUM2=SUM2+ARR(I,K,IX2F)*TMP(K,J)+ARR(I,K,IX2R)*ARR(J,K,IX1F)
101       CONTINUE
          HZERO(IO3+I,IO3+J)=SUM1
          HZERO(IO4+I,IO3+J)=SUM2
102     CONTINUE
100   CONTINUE
      DO 110 I=1,NRROOT
        DO 111 J=1,NRROOT
          TMP(I,J)=ESMALL(I)*ARR(J,I,IX2F)+ARR(J,I,IX2R)
111     CONTINUE
110   CONTINUE
      DO 120 I=1,NRROOT
        DO 122 J=1,NRROOT
          SUM2=HZERO(IO4+I,IO4+J)
          DO 121 K=1,NRROOT
            SUM2=SUM2+ARR(I,K,IX2F)*TMP(K,J)+
     *                ARR(I,K,IX2R)*ARR(J,K,IX2F)
121       CONTINUE
          HZERO(IO4+I,IO4+J)=SUM2
122     CONTINUE
120   CONTINUE
      DO 200 I1=1,NRROOT
        I2=I1+NRROOT
        I3=I1+2*NRROOT
        I4=I1+3*NRROOT
        DO 201 J1=1,NRROOT
          J2=J1+NRROOT
          J3=J1+2*NRROOT
          J4=J1+3*NRROOT
          HZERO(I1,J2)=HZERO(J2,I1)
          HZERO(I1,J3)=HZERO(J3,I1)
          HZERO(I1,J4)=HZERO(J4,I1)
          HZERO(I2,J3)=HZERO(J3,I2)
          HZERO(I2,J4)=HZERO(J4,I2)
          HZERO(I3,J4)=HZERO(J4,I3)
          SZERO(I1,J2)=SZERO(J2,I1)
          SZERO(I1,J3)=SZERO(J3,I1)
          SZERO(I1,J4)=SZERO(J4,I1)
          SZERO(I2,J3)=SZERO(J3,I2)
          SZERO(I2,J4)=SZERO(J4,I2)
          SZERO(I3,J4)=SZERO(J4,I3)
201     CONTINUE
200   CONTINUE
      RETURN
      END
