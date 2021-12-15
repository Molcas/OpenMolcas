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
      SUBROUTINE AIJK(ITAI,L0,L1,L2,L3)
      IMPLICIT REAL*8 (A-H,O-Z)
      External Int8, Int2, int1, int4
#include "SysDef.fh"
#include "files_guga.fh"
      DIMENSION ITAI(*),L0(*),L1(*),L2(*),L3(*)
#include "real_guga.fh"
#include "integ.fh"
      COMMON/CNSTS/D0,D1,D2
#include "addr_guga.fh"
*
      IOUT=0
      NMAT=0
      L=0
      DO 10 NI=1,LN
      DO 20 NJ=1,NI
      I=ICH(NI)
      J=ICH(NJ)
      IF(I.GT.J)GO TO 29
      I=ICH(NJ)
      J=ICH(NI)
29    DO 30 NK=1,LN
      IOUT=IOUT+1
      ICOP1(IOUT)=0
      IF(IOUT.LT.NBUF)GO TO 460
      ICOP1(nCOP+1)=NBUF
      CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT=NMAT+NBUF
      IOUT=0
460   K=ICH(NK)
      IOUT=IOUT+1
*      IND=I+2**10*J
      IND=IOR(I,ISHFT(J,10))
*      ICOP1(IOUT)=IND+2**20*K
      ICOP1(IOUT)=IOR(IND,ISHFT(K,20))
      IF(IOUT.LT.NBUF)GO TO 31
      ICOP1(nCOP+1)=NBUF
      CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT=NMAT+NBUF
      IOUT=0
C     DOUBLET-VALENCE INTERACTIONS
31    ITURN=1
      ITT1=1
      ITT2=0
150   IT1=ITT1*MXVERT
      IT2=ITT2*MXVERT
      JJ=0
      IF(ITT1.NE.0)JJ=IRC(ITT1)
      JJD=0
      IF(ITT1.NE.0)JJD=JRC(ITT1)
      II=0
      IF(ITT2.NE.0)II=IRC(ITT2)
      IID=0
      IF(ITT2.NE.0)IID=JRC(ITT2)
      IF(I.NE.K)GO TO 520
      IF(J.EQ.K)GO TO 524
      CALL INT61(L,J,I,IT1,IT2,II,IID,JJ,JJD,ITURN,ITAI,
     *L0,L1,L2,L3)
      CALL INT62(L,J,I,IT1,IT2,II,IID,JJ,JJD,ITURN,ITAI,
     *L0,L1,L2,L3)
      GO TO 35
524   CALL INT9(L,K,I,IT1,IT2,II,IID,JJ,JJD,ITURN,ITAI,
     *L0,L1,L2,L3)
      GO TO 35
520   IF(J.NE.K)GO TO 535
      CALL INT4(L,K,I,IT1,IT2,II,IID,JJ,JJD,ITURN,ITAI,
     *L0,L1,L2,L3)
      GO TO 35
535   IF(I.NE.J)GO TO 546
      CALL INT8(L,K,I,IT1,IT2,II,IID,JJ,JJD,ITURN,ITAI,
     *L0,L1,L2,L3)
      GO TO 35
546   IF(K.LT.I)GO TO 540
      CALL INT3(J,I,L,K,IT1,IT2,II,IID,JJ,JJD,ITURN,ITAI,
     *L0,L1,L2,L3)
      GO TO 35
540   IF(K.LT.J)GO TO 545
      CALL INT2(L,K,J,I,IT1,IT2,II,IID,JJ,JJD,ITURN,ITAI,
     *L0,L1,L2,L3)
      GO TO 35
545   CALL INT1(L,K,J,I,IT1,IT2,II,IID,JJ,JJD,ITURN,ITAI,
     *L0,L1,L2,L3)
35    GO TO (101,102,103,104,105,30),ITURN
C     TRIPLET-DOUBLET INTERACTIONS
101   ITURN=2
      ITT1=2
      ITT2=1
      GO TO 150
C     SINGLET-DOUBLET INTERACTIONS
102   ITURN=3
      ITT1=3
      ITT2=1
      GO TO 150
C     VALENCE-DOUBLET INTERACTIONS
103   ITURN=4
      ITT1=0
      ITT2=1
      GO TO 150
C     DOUBLET-TRIPLET INTERACTIONS
104   ITURN=5
      ITT1=1
      ITT2=2
      GO TO 150
C     DOUBLET-SINGLET INTERACTIONS
105   ITURN=6
      ITT1=1
      ITT2=3
      GO TO 150
30    CONTINUE
20    CONTINUE
10    CONTINUE
      ICOP1(nCOP+1)=IOUT
      CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT=NMAT+IOUT
      ICOP1(nCOP+1)=-1
      CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      WRITE(IW,600)NMAT
600   FORMAT(/6X,'COEFFICIENTS FOR AIJK',I9)
      RETURN
      END
