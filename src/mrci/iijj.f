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
      SUBROUTINE IIJJ(ICSPCK,INTSYM,HDIAG,FC,FIIJJ,FIJIJ)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
      DIMENSION ICSPCK(*),INTSYM(*),HDIAG(*),
     *          FC(*),FIIJJ(*),FIJIJ(*)
      DIMENSION IOC(55)
*
      JO(L)=ICUNP(ICSPCK,L)
      JSYM(L)=JSUNP(INTSYM,L)
*
      IAD27=0
      II1=0
      ILIM=4
      IF(IFIRST.NE.0)ILIM=2
      IRL=IRC(ILIM)
      DO 100 IR=1,IRL
      DO 110 I=1,LN
      II1=II1+1
      JOJ=JO(II1)
      IF(JOJ.GT.1)JOJ=JOJ-1
      IOC(I)=JOJ
110   CONTINUE
      NSS=MUL(JSYM(IR),LSYM)
      SUM=0.0D00
      DO 111 I=1,LN
        IJ=IROW(I)
        IF(IOC(I).EQ.0)GO TO 111
        DO 113 J=1,I-1
          IJ=IJ+1
          IF(IOC(J).EQ.0)GO TO 113
          SUM=SUM+IOC(I)*(IOC(J)*FIIJJ(IJ)-FIJIJ(IJ))
113     CONTINUE
        IJ=IJ+1
        SUM=SUM+(IOC(I)-1)*FIIJJ(IJ)+IOC(I)*FC(IJ)
111   CONTINUE
      IF(IR.GT.IRC(1))GO TO 120
      HDIAG(IR)=SUM
      IF(IR.NE.IRC(1))GO TO 100
      CALL dDAFILE(Lu_27,1,HDIAG,IRC(1),IAD27)
      GO TO 100
120   IND=0
      IF(IR.GT.IRC(2)) GOTO 130
      NA1=NVIRP(NSS)+1
      NA2=NVIRP(NSS)+NVIR(NSS)
      IF(NA2.LT.NA1)GO TO 100
      DO 121 NA=NA1,NA2
      IND=IND+1
      IA=IROW(LN+NA)
      SUM1=SUM+FC(IA+LN+NA)
      DO 122 I=1,LN
      IF(IOC(I).EQ.0)GO TO 122
      SUM1=SUM1+IOC(I)*FIIJJ(IA+I)-FIJIJ(IA+I)
122   CONTINUE
      HDIAG(IND)=SUM1
121   CONTINUE
      CALL dDAFILE(Lu_27,1,HDIAG,IND,IAD27)
      GO TO 100
130   IND=0
      DO 141 NA=1,NVIRT
      NSA=MUL(NSS,NSM(LN+NA))
      NB1=NVIRP(NSA)+1
      NB2=NVIRP(NSA)+NVIR(NSA)
      IF(NB2.GT.NA)NB2=NA
      IF(NB2.LT.NB1)GO TO 141
      IA=IROW(LN+NA)
      IAV=IA+LN
      DO 142 NB=NB1,NB2
      IND=IND+1
      IB=IROW(LN+NB)
      IBV=IB+LN
      TERM=SUM+FIIJJ(IAV+NB)+FC(IAV+NA)+FC(IBV+NB)
      IF(IR.LE.IRC(3)) THEN
        SUM1=TERM-FIJIJ(IAV+NB)
      ELSE
        SUM1=TERM+FIJIJ(IAV+NB)
      END IF
      DO 143 I=1,LN
      IF(IOC(I).EQ.0)GO TO 143
      TERM=IOC(I)*(FIIJJ(IA+I)+FIIJJ(IB+I))-FIJIJ(IA+I)-FIJIJ(IB+I)
      SUM1=SUM1+TERM
143   CONTINUE
      HDIAG(IND)=SUM1
142   CONTINUE
141   CONTINUE
      IF(IND.GT.0)CALL dDAFILE(Lu_27,1,HDIAG,IND,IAD27)
100   CONTINUE
      RETURN
      END
