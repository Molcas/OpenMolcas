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
      SUBROUTINE ONEEL_GUGA()
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "real_guga.fh"
#include "integ.fh"
#include "files_guga.fh"
      COMMON/CNSTS/D0,D1,D2
#include "addr_guga.fh"
*
      IOUT=0
      NMAT=0
      ITYP=0
      DO 10 NK=1,LN
      DO 20 NI=1,NK
      K=ICH(NK)
      I=ICH(NI)
      IF(K.GT.I)GO TO 19
      K=ICH(NI)
      I=ICH(NK)
19    NSK=NSM(K)
      KJS=IJ(K+1)+1
      KJL=IJ(K)
      NSI=NSM(I)
      IF(NSI.NE.NSK)GO TO 20
      IOUT=IOUT+1
      ICOP1(IOUT)=0
      IF(IOUT.LT.NBUF)GO TO 460
      ICOP1(nCOP+1)=NBUF
      CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT=NMAT+NBUF
      IOUT=0
460   IOUT=IOUT+1
*      ICOP1(IOUT)=I+2**10*K
      ICOP1(IOUT)=IOR(I,ISHFT(K,10))
      IF(IOUT.LT.NBUF)GO TO 21
      ICOP1(nCOP+1)=NBUF
      CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT=NMAT+NBUF
      IOUT=0
21    IF(I.EQ.K)GO TO 20
      DO 101 ITT=1,ILIM
      IT1=(ITT-1)*MXVERT
      IT2=IT1
      DO 30 J=KJS,KJL
      IWAY(K)=1
32    KM=K
      J2(KM+1)=J
      J1(KM+1)=J
      CALL LOOP1(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.1)GO TO 30
41    KM=KM-1
      IWAY(KM)=1
      IF(KM.EQ.I)GO TO 51
42    CALL LOOP5(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.0)GO TO 41
      KM=KM+1
      IF(KM.EQ.K)GO TO 32
      GO TO 42
51    IWAY(I)=1
52    KM=I
      CALL LOOP3(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.1)GO TO 53
      CALL COMP(I,J,ITYP,I,IT1,IT2)
      GO TO 52
53    KM=KM+1
      IF(KM.EQ.K)GO TO 32
      GO TO 42
30    CONTINUE
101   CONTINUE
20    CONTINUE
10    CONTINUE
      ICOP1(nCOP+1)=IOUT
      CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT=NMAT+IOUT
      ICOP1(nCOP+1)=-1
      CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      WRITE(IW,100)NMAT
100   FORMAT(/6X,'COEFFICIENTS FOR IJ',I11)
      RETURN
      END
