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
      SUBROUTINE AI(JTYP,ITAI,L0,L1,L2,L3)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
      DIMENSION ITAI(*),L0(*),L1(*),L2(*),L3(*)
#include "real_guga.fh"
#include "integ.fh"
#include "files_guga.fh"
      COMMON/CNSTS/D0,D1,D2
#include "addr_guga.fh"
      COMMON/D/JNDX(500 000)
      IOUT=0
      NMAT=0
      JMAX=0
      DO 10 NI=1,LN
      IOUT=IOUT+1
      ICOP1(IOUT)=0
      JOUT=0
      IF(IOUT.LT.NBUF)GO TO 460
      ICOP1(NCOP+1)=NBUF
      CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT=NMAT+NBUF
      IOUT=0
460   I=ICH(NI)
      IOUT=IOUT+1
      ICOP1(IOUT)=I
      IF(IOUT.LT.NBUF)GO TO 11
      ICOP1(NCOP+1)=NBUF
      CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT=NMAT+NBUF
      IOUT=0
11    IJS=IJ(I+1)+1
      IJM=IJ(I)
      IF(JTYP.EQ.1)GO TO 101
C     DOUBLET-VALENCE INTERACTIONS
      ITURN=1
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
      ITYP=ITURN
      DO 30 IJJ=IJS,IJM
      ITAIL=IX(IT2+IJJ)
      CALL TAIL(I,IJJ,ITAI,ITAIL,L0,L1,L2,L3,IT1,IT2)
      IWAY(I)=1
32    KM=I
      J2(KM+1)=IJJ
      J1(KM+1)=IJJ
      CALL LOOP1(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.1)GO TO 30
41    KM=KM-1
      IF(KM.EQ.0)GO TO 51
      IWAY(KM)=1
42    CALL LOOP5(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.0)GO TO 41
53    KM=KM+1
      IF(KM.EQ.I)GO TO 32
      GO TO 42
51    DO 80 IN=1,ITAIL
      ICP1=ICOUP(1)+IN
      JND1=JNDX(II+ICP1)
      IF(JND1.EQ.0)GO TO 80
      ICP1=JND1-IID
      IN2=ITAI(IN)
      IF(IN2.EQ.0)GO TO 80
      ICP2=ICOUP1(1)+IN2
      JND2=JNDX(JJ+ICP2)
      IF(JND2.EQ.0)GO TO 80
      ICP2=JND2-JJD
      IOUT=IOUT+1
      JOUT=JOUT+1
      IF(JOUT.GT.JMAX)JMAX=JOUT
      COP(IOUT)=COUP(1)
CPAM96      IND=IOR(ITYP,ISHFT(ICP2,6))
*      IND=ITYP+2**6*ICP2
      IND=IOR(ITYP,ISHFT(ICP2,6))
CPAM96      ICOP1(IOUT)=IOR(IND,ISHFT(ICP1,19))
*      ICOP1(IOUT)=IND+2**19*ICP1
      ICOP1(IOUT)=IOR(IND,ISHFT(ICP1,19))
      IF(IOUT.LT.NBUF)GO TO 80
      ICOP1(NCOP+1)=NBUF
      CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT=NMAT+NBUF
      IOUT=0
80    CONTINUE
      GO TO 53
30    CONTINUE
      GO TO (101,102,10),ITURN
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
10    CONTINUE
      ICOP1(NCOP+1)=IOUT
      CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT=NMAT+IOUT
      ICOP1(NCOP+1)=-1
      CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      CHKSUM=0.0D0
      DO I=1,NCOP
       CHKSUM=CHKSUM+COP(I)
      END DO
      CALL ADD_INFO('GUGA_CHKSUM',[CHKSUM],1,8)
      IF(JTYP.EQ.0)WRITE(IW,600)NMAT
600   FORMAT(/,6X,'COEFFICIENTS FOR AI',I11)
      IF(JTYP.EQ.0)THEN
        RETURN
      END IF
      WRITE(IW,601)NMAT
601   FORMAT(/,6X,'COEFFICIENTS FOR ABCI',I9)
      IAD10(1)=JMAX
      WRITE(IW,602)JMAX
602   FORMAT(6X,'MAXIMUM NUMBER OF ELEMENTS',I6)
      RETURN
      END
