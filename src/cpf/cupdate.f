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
      SUBROUTINE CUPDATE(JSY,INDEX,C,S,AP,BST,T2,ENP)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION JSY(*),INDEX(*),C(*),S(*),AP(*),BST(*),
     *          T2(*),ENP(*)
#include "SysDef.fh"
#include "cpfmcpf.fh"
#include "files_cpf.fh"
*
      JSYM(L)=JSUNP(JSY,L)
C

      W=WLEV

C VALENCE
      IP=IRC(1)
      DO 6 I=1,IP
        C(I)=AP(I)*C(I)-S(I)
6     CONTINUE

C SINGLES
      IP=IRC(2)-IRC(1)
      IN=IRC(1)
      DO 5 I=1,IP
        NS1=JSYM(IN+I)
        NSIL=MUL(NS1,LSYM)
        INUM=NVIR(NSIL)
        IST=INDEX(IN+I)+1
        CALL VSMSB(C(IST),1,AP(IN+I),S(IST),1,C(IST),1,INUM)
5     CONTINUE

C DOUBLES
      IP=IRC(4)-IRC(2)
      IN=IRC(2)
      DO 10 I=1,IP
        NS1=JSYM(IN+I)
        NSIL=MUL(NS1,LSYM)
        INUM=NNS(NSIL)
        IST=INDEX(IN+I)+1
        CALL VSMSB(C(IST),1,AP(IN+I),S(IST),1,C(IST),1,INUM)
10    CONTINUE

C WRITE GRADIENT ONTO DISK, UNIT=30
      IAD=IADDP(ITPUL)
      CALL dDAFILE(Lu_30,1,C,NCONF,IAD)

C Reuse array S for HCOUT that was written in IJIJ.
      IAD=IAD25S
      DO 77 III=1,NCONF,nCOP
        JJJ=MIN(nCOP,NCONF+1-III)
        CALL dDAFILE(Lu_25,2,S(III),JJJ,IAD)
77    CONTINUE

C VALENCE
      IP=IRC(1)
      DO 7 I=1,IP
        APW=W-AP(I)
        T2(I)=S(I)+APW
        C(I)=C(I)/T2(I)
        EMP=SQRT(ENP(I))
        C(I)=C(I)*EMP
        IF(I.EQ.IREF0)C(I)=0.0D0
7     CONTINUE

C SINGLES
      IP=IRC(2)-IRC(1)
      IN=IRC(1)
      DO 8 I=1,IP
        NS1=JSYM(IN+I)
        NSIL=MUL(NS1,LSYM)
        INUM=NVIR(NSIL)
        IST=INDEX(IN+I)+1
        APW=W-AP(IN+I)
        CALL VSADD(S(IST),1,APW,T2,1,INUM)
        CALL VDIV(T2,1,C(IST),1,C(IST),1,INUM)
        EMP=SQRT(ENP(IN+I))
        CALL VSMUL(C(IST),1,EMP,C(IST),1,INUM)
8     CONTINUE

C DOUBLES
      IP=IRC(4)-IRC(2)
      IN=IRC(2)
      DO 11 I=1,IP
        NS1=JSYM(IN+I)
        NSIL=MUL(NS1,LSYM)
        INUM=NNS(NSIL)
        IST=INDEX(IN+I)+1
        APW=W-AP(IN+I)
        CALL VSADD(S(IST),1,APW,T2,1,INUM)
        CALL VDIV(T2,1,C(IST),1,C(IST),1,INUM)
        EMP=SQRT(ENP(IN+I))
        CALL VSMUL(C(IST),1,EMP,C(IST),1,INUM)
11    CONTINUE
C
      IAD=IADDP(ITPUL+1)
      CALL dDAFILE(Lu_CI,1,C,NCONF,IAD)
      IADDP(ITPUL+2)=IAD
      IF(IPRINT.GE.15)WRITE(6,999)(C(I),I=1,NCONF)
999   FORMAT(6X,'C(UPD)',5F10.6)
      A=DDOT_(NCONF,C,1,C,1)
      IF(A.GT.D2) THEN
        WRITE(6,*)'CUPDATE Error: A>2.0D0 (See code.)'
        CALL Abend
      END IF
      IF(ITPUL.EQ.1)BST(1)=A
      RETURN
      END
