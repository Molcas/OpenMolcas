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
      SUBROUTINE EMPTY(BUF,IBUF,LASTAD,SO,KBUF,NTPB)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "SysDef.fh"
#include "real_guga.fh"
#include "integ.fh"
#include "files_guga.fh"
      DIMENSION BUF(kBuf),IBUF(kBuf+2),LASTAD(*),SO(*)
      COMMON/CNSTS/D0,D1,D2
#include "addr_guga.fh"
      COMMON/D/JNDX(500 000)
*
      JO(L)=ICUNP(ICASE,L)
*
*
      ISUM=JRC(ILIM)
      IOUT=0
      NMAT=0
      ICLR=NTPB
      IN=ICLR+1
      NBX=0
      IOFF=0
      DO 10 II=1,ISUM
      IF(II.GT.JRC(1))GO TO 11
      IVL=IV0
      KK=II
      GO TO 15
11    IF(II.GT.JRC(2))GO TO 12
      IVL=IV1
      KK=II-JRC(1)
      GO TO 15
12    IF(II.GT.JRC(3))GO TO 13
      IVL=IV2
      KK=II-JRC(2)
      GO TO 15
13    IVL=IV3
      KK=II-JRC(3)
15    IOUT=IOUT+1
      ICOP1(IOUT)=0
      IF(IOUT.LT.NBUF)GO TO 460
      ICOP1(nCOP+1)=NBUF
      CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT=NMAT+NBUF
      IOUT=0
460   IVL0=IV0-IVL
*      IND=KK+2**16*IVL0
      IND=IOR(KK,ISHFT(IVL0,16))
      IOUT=IOUT+1
      ICOP1(IOUT)=IND
      IF(IOUT.LT.NBUF)GO TO 16
      ICOP1(nCOP+1)=NBUF
      CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT=NMAT+NBUF
      IOUT=0
16    IJJ=0
      JJ=(II-1)*LN
      DO 20 I=1,LN
      JJ1=JO(JJ+I)
      DO 25 J=1,I
      IJJ=IJJ+1
      IN=IN+1
      IF(IN.LE.ICLR)GO TO 100
      IN=1
      DO 105 IIQQ=1,ICLR
      SO(IIQQ)=D0
105   CONTINUE
      NBX=NBX+1
      IADR=LASTAD(NBX)
110   CONTINUE
      IF(IADR.EQ.-1)GO TO 120
C     FPS
      CALL dDAFILE(Lu_11,2, BUF,KBUF,  IADR)
      CALL iDAFILE(Lu_11,2,IBUF,KBUF+2,IADR)
      LENGTH=IBUF(KBUF+1)
      IADR  =IBUF(KBUF+2)
      IF(LENGTH.EQ.0)GO TO 110
      DO 130 IIQQ=1,LENGTH
      IQ=IBUF(IIQQ)-IOFF
      SO(IQ)=BUF(IIQQ)
130   CONTINUE
      GO TO 110
120   IOFF=IOFF+ICLR
100   IF(JJ1.EQ.0)GO TO 25
      JJ2=JO(JJ+J)
      IF(JJ2.EQ.0)GO TO 25
      ITYP=0
      IF(I.EQ.J)ITYP=1
      IKK=IJJ
      IF(ITYP.EQ.1)IKK=I
      IF(SO(IN).EQ.D0)GO TO 25
      IOUT=IOUT+1
CPAM96      ICOP1(IOUT)=IOR(ITYP,ISHFT(IKK,1))
*      ICOP1(IOUT)=ITYP+2*IKK
      ICOP1(IOUT)=IOR(ITYP,ISHFT(IKK,1))
      COP(IOUT)=SO(IN)
      IF(IOUT.LT.NBUF)GO TO 25
      ICOP1(nCOP+1)=NBUF
      CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT=NMAT+NBUF
      IOUT=0
25    CONTINUE
20    CONTINUE
10    CONTINUE
      RETURN
      END
