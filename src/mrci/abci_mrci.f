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
      SUBROUTINE ABCI_MRCI(INTSYM,indx,C,S,BMN,IBMN,BIAC,BICA,BUFIN)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
      DIMENSION INTSYM(*),indx(*),C(*),S(*),BMN(*),IBMN(*),
     *BIAC(*),BICA(*),BUFIN(*)
*
      JSYM(L)=JSUNP(INTSYM,L)
*
      CALL CSCALE(indx,INTSYM,C,SQ2)
      CALL CSCALE(indx,INTSYM,S,SQ2INV)
      ICHK=0
      INSIN=KBUFF1
      IAD15=IADABCI
      IADD10=IAD10(4)
      CALL dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
      CALL iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
      LEN=ICOP1(nCOP+1)
      IN=2
      NSAVE=ICOP1(IN)
100   NI=NSAVE
      IOUT=0
110   IN=IN+1
      IF(IN.LE.LEN)GO TO 15
      CALL dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
      CALL iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
      LEN=ICOP1(nCOP+1)
      IF(LEN.LE.0)GO TO 5
      IN=1
15    IF(ICHK.NE.0)GO TO 460
      IF(ICOP1(IN).EQ.0)GO TO 10
      IOUT=IOUT+1
      BMN(IOUT)=COP(IN)
      IBMN(IOUT)=ICOP1(IN)
      GO TO 110
10    ICHK=1
      GO TO 110
460   ICHK=0
      NSAVE=ICOP1(IN)
5     CONTINUE
      DO 20 NB=1,NVIRT
      NSIB=MUL(NSM(LN+NB),NSM(NI))
      NSLB=MUL(NSM(LN+NB),LSYM)
      LB=NB-NVIRP(NSM(LN+NB))
      INS=NVPAIR(NSIB)
      ILOOP=0
72    INSB=INS
73    IF(INSIN.LT.KBUFF1)GO TO 75
      CALL dDAFILE(Lu_70,2,BUFIN,KBUFF1,IAD15)
      INSIN=0
75    INB=KBUFF1-INSIN
      INUMB=INSB
      IF(INSB.GT.INB)INUMB=INB
      IST=INS-INSB+1
      IF(ILOOP.EQ.0)CALL DCOPY_(INUMB,BUFIN(INSIN+1),1,BIAC(IST),1)
      IF(ILOOP.EQ.1)CALL DCOPY_(INUMB,BUFIN(INSIN+1),1,BICA(IST),1)
      INSIN=INSIN+INUMB
      INSB=INSB-INUMB
      IF(INSB.GT.0)GO TO 73
      ILOOP=ILOOP+1
      IF(ILOOP.EQ.1)GO TO 72
      DO 25 IT=1,IOUT
      IND=IBMN(IT)
*      ICP1=MOD(IND/2**19,2**13)
      ICP1=IBITS(IND,19,13)
      INDA=IRC(1)+ICP1
      IF(JSYM(INDA).NE.NSLB)GO TO 25
      MA=indx(INDA)+LB
*      ICP2=MOD(IND/2**6,2**13)
*      ITYP=MOD(IND,2**6)
      ICP2=IBITS(IND,6,13)
      ITYP=IBITS(IND,0,6)
      IF(INS.EQ.0)GO TO 25
      COPL=BMN(IT)*C(MA)
      INDB=IRC(ITYP)+ICP2
      ICCB=indx(INDB)+1
      IF(ITYP.EQ.3)GO TO 26
      TERM=DDOT_(INS,C(ICCB),1,BICA,1)
      CALL DAXPY_(INS,COPL,BICA,1,S(ICCB),1)
      GO TO 27
26    TERM=DDOT_(INS,C(ICCB),1,BIAC,1)
      CALL DAXPY_(INS,COPL,BIAC,1,S(ICCB),1)
27    S(MA)=S(MA)+BMN(IT)*TERM
25    CONTINUE
20    CONTINUE
      IF(LEN.GE.0)GO TO 100
      CALL CSCALE(indx,INTSYM,C,SQ2INV)
      CALL CSCALE(indx,INTSYM,S,SQ2)
      RETURN
      END
