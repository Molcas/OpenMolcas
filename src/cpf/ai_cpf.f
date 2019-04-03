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
*               1986, Margareta R. A. Blomberg                         *
************************************************************************
      SUBROUTINE AI_CPF(JSY,INDEX,C,S,FC,BUFIN,IBUFIN,A,B,FK,DBK,
     *ENP,EPP,KTYP)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "cpfmcpf.fh"
#include "files_cpf.fh"
      DIMENSION JSY(*),INDEX(*),C(*),S(*),FC(*),BUFIN(*),IBUFIN(*),
     &           A(*),B(*),FK(*),DBK(*),ENP(*),EPP(*)
      PARAMETER (IPOW10=2**10, IPOW20=2**20)
      PARAMETER (IPOW6=2**6, IPOW19=2**19)
*
C     KTYP=0  ,  (A/I)   INTEGRALS
C     KTYP=1  ,  (AI/JK) INTEGRALS
      DIMENSION IPOB(9)
*
      JSYM(L)=JSUNP_CPF(JSY,L)
*
      NK = 0 ! dummy initialize
      NSK= 0 ! dummy initialize
      INUM=IRC(4)-IRC(3)
      CALL PSQ2(C,S,MUL,INDEX,JSY,NDIAG,INUM,IRC(3),LSYM,NVIRT,SQ2)
      NVT=IROW(NVIRT+1)
      ICHK=0
      IJOLD=0
      NOB2=IROW(NORBT+1)
      NOT2=IROW(LN+1)
      NOTT=2*NOT2
      NOVST=LN*NVIRT+1+NVT
      LBUF0=RTOI*LBUF
      LBUF1=LBUF0+LBUF+1
      LBUF2=LBUF1+1
      IF(KTYP.EQ.0)IADD10=IAD10(9)
      IF(KTYP.EQ.1)IADD10=IAD10(7)
100   CALL dDAFILE(Lu_CIGuga,2,COP,nCOP,IADD10)
      CALL iDAFILE(Lu_CIGuga,2,iCOP1,nCOP+1,IADD10)
      LEN=ICOP1(nCOP+1)
      IF(LEN.EQ.0)GO TO 100
      IF(LEN.LT.0)GO TO 200
      DO 10 II=1,LEN
      IND=ICOP1(II)
      IF(ICHK.NE.0)GO TO 460
      IF(IND.NE.0)GO TO 11
      ICHK=1
      GO TO 10
460   ICHK=0
      ITURN=0
      IF(IDENS.EQ.1.AND.IJOLD.NE.0)GO TO 20
21    ITURN=1
      IF(KTYP.EQ.1)GO TO 9
      NK=IND
      IJOLD=NK
      NSK=NSM(NK)
      GO TO 20
9     INDI=IND
*      NI=MOD(INDI,1024)
*      NJ=MOD(INDI/IPOW10,1024)
*      NK=MOD(INDI/IPOW20,1024)
      NI=IBITS(INDI,0,10)
      NJ=IBITS(INDI,10,10)
      NK=IBITS(INDI,20,10)
      NSIJ=MUL(NSM(NI),NSM(NJ))
      NSK=MUL(NSIJ,NSM(NK))
      IJ=IROW(NI)+NJ
      IF(IJ.EQ.IJOLD)GO TO 20
      IJOLD=IJ
      IADR=LASTAD(NOVST+NOTT+IJ)
      DO 105 INN=1,NOB2
      FC(INN)=D0
105   CONTINUE
90    CALL iDAFILE(Lu_TiABIJ,2,IBUFIN,LBUF2,IADR)
      LENGTH=IBUFIN(LBUF1)
      IADR=IBUFIN(LBUF2)
      IF(LENGTH.EQ.0)GO TO 91
      CALL SCATTER(LENGTH,FC,IBUFIN(LBUF0+1),BUFIN)
91    IF(IADR.NE.-1) GO TO 90
C     FORM VECTOR FK
20    NA1=NSYS(NSK)+1
      NA2=NSYS(NSK+1)
      INK=0
      IF(NA2.LT.NA1)GO TO 10
      DO 13 NA=NA1,NA2
      INK=INK+1
      NAK=IROW(LN+NA)+NK
      IF(ITURN.EQ.0)FC(NAK)=FK(INK)
      IF(ITURN.EQ.1)FK(INK)=FC(NAK)
13    CONTINUE
      IF(ITURN.EQ.0)GO TO 21
      GO TO 10
11    IF(INK.EQ.0)GO TO 10
CPAM97      ITYP=IAND(IND,63)
CPAM97      ICP2=IAND(ISHFT(IND,-6),8191)
CPAM97      ICP1=IAND(ISHFT(IND,-19),8191)
*      ITYP=MOD(IND,64)
*      ICP2=MOD(IND/IPOW6,8192)
*      ICP1=MOD(IND/IPOW19,8192)
      ITYP=IBITS(IND, 0, 6)
      ICP2=IBITS(IND,6,13  )
      ICP1=IBITS(IND,19,13  )
      IF(ITYP.GT.1)GO TO 12
      INDA=ICP1
      INDB=IRC(1)+ICP2
      INNY=INDEX(INDB)+1
      IF(IDENS.EQ.1)GO TO 41
      IF(INDA.NE.IREF0)GO TO 42
      COPI=COP(II)/SQRT(ENP(INDB))
      CALL DAXPY_(INK,COPI,FK,1,S(INNY),1)
      IF(ITER.EQ.1)GO TO 10
      TERM=DDOT_(INK,FK,1,C(INNY),1)
      EPP(INDB)=EPP(INDB)+COPI*TERM
      GO TO 10
42    FACS=D1
      COPI=COP(II)*C(INDA)
      CALL DAXPY_(INK,COPI*FACS,FK,1,S(INNY),1)
      TERM=DDOT_(INK,FK,1,C(INNY),1)
      S(INDA)=S(INDA)+COP(II)*FACS*TERM
      GO TO 10
41    IF(INDA.EQ.IREF0)COPI=C(INDA)*COP(II)/ENP(INDB)
      IF(INDA.NE.IREF0)COPI=C(INDA)*COP(II)/
     *(SQRT(ENP(INDA))*SQRT(ENP(INDB)))
      CALL DAXPY_(INK,COPI,C(INNY),1,FK,1)
C      WRITE(6,654)NK,NSK,INDB
C  654 FORMAT(1X,'TYP1,NK,NSK,INDB',3I7)
C      WRITE(6,653)(FK(I),I=1,INK)
C  653 FORMAT(1X,'FK',5F12.6)
      GO TO 10
12    IF(ITER.EQ.1)GO TO 10
      INDA=IRC(1)+ICP1
      INDB=IRC(ITYP)+ICP2
      INMY=INDEX(INDA)+1
      INNY=INDEX(INDB)+1
      MYSYM=JSYM(INDA)
      NYSYM=MUL(MYSYM,NSK)
      MYL=MUL(MYSYM,LSYM)
      NYL=MUL(NYSYM,LSYM)
      IFT=0
      IF(ITYP.EQ.2)IFT=1
      CALL IPO_CPF(IPOB,NVIR,MUL,NSYM,NYL,IFT)
      NVM=NVIR(MYL)
      IF(IDENS.EQ.1)GO TO 210
      FACS=D1
      CALL SETZ(DBK,INK)
      CALL DAXPY_(INK,COP(II),FK,1,DBK,1)
      IF(NYL.NE.1)GO TO 25
      IF(IFT.EQ.0)CALL SQUAR_CPF(C(INNY+IPOB(MYL)),A,NVM)
      IF(IFT.EQ.1)CALL SQUARM_CPF(C(INNY+IPOB(MYL)),A,NVM)
      CALL SETZ(B,NVM)
      CALL FMMM(DBK,A,B,1,NVM,INK)
      CALL DAXPY_(NVM,FACS,B,1,S(INMY),1)
      SIGN=D1
      IF(IFT.EQ.1)SIGN=-D1
      IOUT=INNY+IPOB(MYL)-1
      DO 125 I=1,NVM
      DO 130 J=1,I
      IOUT=IOUT+1
      TERM=DBK(I)*C(INMY+J-1)+SIGN*DBK(J)*C(INMY+I-1)
      S(IOUT)=S(IOUT)+FACS*TERM
130   CONTINUE
      IF(IFT.EQ.1)GO TO 125
      TERM=DBK(I)*C(INMY+I-1)
      S(IOUT)=S(IOUT)-FACS*TERM
125   CONTINUE
      GO TO 10
25    NKM=INK*NVM
      CALL SETZ(B,NVM)
      IF(NSK.GT.MYL)GO TO 26
      IF(IFT.EQ.1)CALL VNEG_CPF(DBK,1,DBK,1,INK)
      CALL FMMM(DBK,C(INNY+IPOB(MYL)),B,1,NVM,INK)
      CALL DAXPY_(NVM,FACS,B,1,S(INMY),1)
      CALL SETZ(B,NKM)
      CALL FMMM(DBK,C(INMY),B,INK,NVM,1)
      CALL DAXPY_(NKM,FACS,B,1,S(INNY+IPOB(MYL)),1)
      GO TO 10
26    CALL FMMM(C(INNY+IPOB(NSK)),DBK,B,NVM,1,INK)
      CALL DAXPY_(NVM,FACS,B,1,S(INMY),1)
      CALL SETZ(B,NKM)
      CALL FMMM(C(INMY),DBK,B,NVM,INK,1)
      CALL DAXPY_(NKM,FACS,B,1,S(INNY+IPOB(NSK)),1)
      GO TO 10
210   CALL SETZ(B,INK)
      COPI=COP(II)/(SQRT(ENP(INDA))*SQRT(ENP(INDB)))
C      WRITE(6,652)IFT,NYL,NSK,MYL,INDA,INDB
C  652 FORMAT(1X,'TYP2',6I7)
      IF(NYL.NE.1)GO TO 225
      IF(IFT.EQ.0)CALL SQUAR_CPF(C(INNY+IPOB(MYL)),A,NVM)
      IF(IFT.EQ.1)CALL SQUARN_CPF(C(INNY+IPOB(MYL)),A,NVM)
      CALL FMMM(C(INMY),A,B,1,INK,NVM)
227   CALL VSMA(B,1,COPI,FK,1,FK,1,INK)
C      WRITE(6,651)(FK(I),I=1,INK)
C  651 FORMAT(1X,'FK',5F12.6)
      GO TO 10
225   IF(NSK.GT.MYL)GO TO 226
      CALL FMMM(C(INNY+IPOB(MYL)),C(INMY),B,INK,1,NVM)
      IF(IFT.EQ.1)COPI=-COPI
      GO TO 227
226   CALL FMMM(C(INMY),C(INNY+IPOB(NSK)),B,1,INK,NVM)
      GO TO 227
10    CONTINUE
      GO TO 100
200   IF(IDENS.EQ.0)GO TO 201
      NA1=NSYS(NSK)+1
      NA2=NSYS(NSK+1)
      INK=0
      IF(NA2.LT.NA1)GO TO 201
      DO 213 NA=NA1,NA2
      INK=INK+1
      NAK=IROW(LN+NA)+NK
      FC(NAK)=FK(INK)
213   CONTINUE
201   CALL DSQ2(C,S,MUL,INDEX,JSY,NDIAG,INUM,IRC(3),
     *LSYM,NVIRT,SQ2)
      RETURN
      END
