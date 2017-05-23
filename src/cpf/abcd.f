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
      SUBROUTINE ABCD(JSY,INDEX,ISAB,C,S,ACBDS,ACBDT,BUFIN)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "cpfmcpf.fh"
#include "files_cpf.fh"
      DIMENSION JSY(*),INDEX(*),ISAB(*),C(*),S(*),ACBDS(*),ACBDT(*),
     &           BUFIN(*)
      IAD16=0
      KBUFF1=2*9600
      INSIN=KBUFF1
      INUM=IRC(4)-IRC(3)
      CALL PSQ2(C,S,MUL,INDEX,JSY,NDIAG,INUM,IRC(3),LSYM,NVIRT,SQ2)
      NVT=IROW(NVIRT+1)
      NOV=(NVT-1)/IPASS+1
      IACMAX=0
      DO 70 ISTEP=1,IPASS
      IACMIN=IACMAX+1
      IACMAX=IACMAX+NOV
      IF(IACMAX.GT.NVT)IACMAX=NVT
      IF(IACMIN.GT.IACMAX)GO TO 70
      DO 40 ISYM=1,NSYM
      IST1=IRC(3)+JJS(ISYM+9)+1
      IFIN1=IRC(3)+JJS(ISYM+10)
      INPS=IFIN1-IST1+1
      IST2=IRC(2)+JJS(ISYM)+1
      IFIN2=IRC(2)+JJS(ISYM+1)
      INPT=IFIN2-IST2+1
      ITAIL=INPS+INPT
      IF(ITAIL.EQ.0)GO TO 40
      IN1=-NVIRT
      DO 50 NA=1,NVIRT
      IN1=IN1+NVIRT
      DO 60 NC=1,NA
      IAC=IROW(NA)+NC
      IF(IAC.LT.IACMIN)GO TO 60
      IF(IAC.GT.IACMAX)GO TO 60
      IF(NA.EQ.1)GO TO 60
      NSAC=MUL(NSM(LN+NA),NSM(LN+NC))
      NSACL=MUL(NSAC,LSYM)
      IF(NSACL.NE.ISYM)GO TO 60
      ISAC=ISAB(IN1+NC)
      NDMAX=NSYS(NSM(LN+NC)+1)
      IF(NDMAX.GT.NA)NDMAX=NA
      INS=ISAB(IN1+NDMAX)
      ILOOP=0
72    INSB=INS
73    IF(INSIN.LT.KBUFF1)GO TO 75
      CALL dDAFILE(Lu_TiABCD,2,BUFIN,KBUFF1,IAD16)
      INSIN=0
75    INB=KBUFF1-INSIN
      INUMB=INSB
      IF(INSB.GT.INB)INUMB=INB
      IST=INS-INSB+1
      IF(ILOOP.EQ.0)CALL DCOPY_(INUMB,BUFIN(INSIN+1),1,
     *ACBDS(IST),1)
      IF(ILOOP.EQ.1)CALL DCOPY_(INUMB,BUFIN(INSIN+1),1,
     *ACBDT(IST),1)
      INSIN=INSIN+INUMB
      INSB=INSB-INUMB
      IF(INSB.GT.0)GO TO 73
      ILOOP=ILOOP+1
      IF(ILOOP.EQ.1)GO TO 72
      IF(INPS.EQ.0)GO TO 11
      DO 10 INDA=IST1,IFIN1
      FACS=D1
      TERM=DDOT_(INS,C(INDEX(INDA)+1),1,ACBDS,1)
      S(INDEX(INDA)+ISAC)=S(INDEX(INDA)+ISAC)+FACS*TERM
      CALL DAXPY_(INS,FACS*C(INDEX(INDA)+ISAC),ACBDS,1,
     &                                        S(INDEX(INDA)+1),1)
10    CONTINUE
11    IF(INPT.EQ.0.OR.NA.EQ.NC)GO TO 60
      DO 30 INDA=IST2,IFIN2
      FACS=D1
      TERM=DDOT_(INS,C(INDEX(INDA)+1),1,ACBDT,1)
      S(INDEX(INDA)+ISAC)=S(INDEX(INDA)+ISAC)+FACS*TERM
      CALL DAXPY_(INS,FACS*C(INDEX(INDA)+ISAC),ACBDT,1,
     &                                        S(INDEX(INDA)+1),1)
30    CONTINUE
60    CONTINUE
50    CONTINUE
40    CONTINUE
70    CONTINUE
      CALL DSQ2(C,S,MUL,INDEX,JSY,NDIAG,INUM,IRC(3),LSYM,NVIRT,SQ2)
      RETURN
      END
