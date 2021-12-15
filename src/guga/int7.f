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
* Copyright (C) 1986, Per E. M. Siegbahn                               *
************************************************************************
      SUBROUTINE INT7(I,K,L,IDIAG,BUFOUT,INDOUT,ICAD,IBUFL,
     *KBUF,NTPB)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
      DIMENSION BUFOUT(*),INDOUT(*),ICAD(*),IBUFL(*)
C     I.LT.L  I.EQ.K  J.EQ.L
#include "real_guga.fh"
#include "integ.fh"
#include "files_guga.fh"
      COMMON/CNSTS/D0,D1,D2
#include "addr_guga.fh"
      COMMON/D/JNDX(500 000)
*
      IJJ=0 ! dummy initialize
      KBUF0=RTOI*KBUF
      KBUF1=KBUF0+KBUF+1
      KBUF2=KBUF1+1
      IDIV=RTOI
      ITYP=0
      IF(IDIAG.EQ.1)IJJ=L*(L-1)/2+K
      LJS=IJ(L+1)+1
      LJM=IJ(L)
      DO 100 ITT=1,ILIM
      IT1=(ITT-1)*MXVERT
      IT2=IT1
      DO 10 LJ=LJS,LJM
      ITURN=0
      IF(IDIAG.EQ.1)ITURN=1
33    IWAY(L)=1
32    KM=L
      J2(KM+1)=LJ
      J1(KM+1)=LJ
      JM(KM)=IVF0+1
      JM1(KM)=IVF0+1
      IF(ITURN.EQ.0)CALL LOOP7(KM,ISTOP,IT1,IT2)
      IF(ITURN.EQ.1)CALL LOOP8(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.1)GO TO 54
      IF(IDIAG.EQ.1.AND.J1(KM).NE.J2(KM))GO TO 32
      GO TO 53
54    IF(ITURN.EQ.1)GO TO 10
      ITURN=1
      GO TO 33
53    KM=KM-1
      IWAY(KM)=1
      IF(KM.EQ.I)GO TO 71
62    JM(KM)=IVF0+1
      JM1(KM)=IVF0+1
      IF(ITURN.EQ.0)CALL LOOP17(KM,ISTOP,IT1,IT2)
      IF(ITURN.EQ.1)CALL LOOP21(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.1)GO TO 55
      IF(IDIAG.EQ.1.AND.J1(KM).NE.J2(KM))GO TO 62
      GO TO 53
55    KM=KM+1
      IF(KM.EQ.L)GO TO 32
      GO TO 62
71    KM=I
      IF(ITURN.EQ.0)CALL LOOP14(KM,ISTOP,IT1,IT2)
      IF(ITURN.EQ.1)CALL LOOP18(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.1)GO TO 73
      IF(ABS(COUP(I)).LT.1.D-06)GO TO 71
      IF(IDIAG.EQ.0)GO TO 105
      IF(ICOUP1(I).NE.ICOUP(I))GO TO 71
      GO TO 106
105   IF(ITURN.EQ.0.OR.IWAY(L).EQ.5)GO TO 106
      IF(ICOUP1(I).LT.ICOUP(I))GO TO 71
      IF(ICOUP1(I).EQ.ICOUP(I))GO TO 71
106   IF(ITURN.EQ.0)COUP(I)=COUP(I)/D2
      IF(IDIAG.EQ.1)GO TO 25
      CALL COMP(I,LJ,ITYP,I,IT1,IT2)
      GO TO 71
25    KM=KM-1
      IF(KM.EQ.0)GO TO 26
      IWAY(KM)=1
27    CALL PATH(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.0)GO TO 25
      KM=KM+1
      IF(KM.EQ.I)GO TO 71
      GO TO 27
26    IVL=J2(1)
      ITAIL=IX(IT1+LJ)
      ISUM=IV0-IVL
      ISU=0
      IF(ISUM.NE.0)ISU=IRC(ISUM)
      DO 104 IN=1,ITAIL
        JND1=JNDX(ISU+ICOUP(1)+IN)
        IF(JND1.EQ.0)GO TO 104
        IPOS=(JND1-1)*LNP+IJJ
        NBN=(IPOS-1)/NTPB+1
        IBUFL(NBN)=IBUFL(NBN)+1
        ICQ=ICAD(NBN)
        ICP=ICQ/IDIV+IBUFL(NBN)
        BUFOUT(ICP)=COUP(I)
        ICPP=ICQ+KBUF0+IBUFL(NBN)
        INDOUT(ICPP)=IPOS
        IF(IBUFL(NBN).LT.KBUF)GO TO 104
        INDOUT(ICQ+KBUF1)=KBUF
        IAD110=IADD11
        CALL iDAFILE(Lu_11,1,INDOUT(ICQ+1),KBUF2,IADD11)
        INDOUT(ICQ+KBUF2)=IAD110
        IBUFL(NBN)=0
104   CONTINUE
      IF(I.EQ.1)GO TO 71
      KM=1
      GO TO 27
73    KM=KM+1
      IF(KM.EQ.L)GO TO 32
      GO TO 62
10    CONTINUE
100   CONTINUE
      RETURN
      END
