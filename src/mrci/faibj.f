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
      SUBROUTINE FAIBJ(INTSYM,INDX,C,S,ABIJ,AIBJ,AJBI,BUF,
     *      IBUF,A,B,F,FSEC)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
#include "WrkSpc.fh"
      DIMENSION INTSYM(*),INDX(*),C(*),S(*),
     *          ABIJ(*),AIBJ(*),AJBI(*),
     *          BUF(NBITM3),IBUF(NBITM3+2),
     *          A(*),B(*),F(*),FSEC(*)
      DIMENSION IPOF(9),IPOA(9),IPOB(9)
      External JSUNP
*------
c
cvv this code is a real compiler killer!
c
* POW: Unnecessary but warningstopping initializations
      iTyp=-1234567
      iCoup=-1234567
      iCoup1=-1234567
c      call getmem('test','chec','real',ldum,ndum)
*------
      CALL CSCALE(INDX,INTSYM,C,SQ2)
      CALL CSCALE(INDX,INTSYM,S,SQ2INV)
      ICHK=0
      IFAB=0
      NOVST=LN*NVIRT+1+(NVIRT*(NVIRT+1))/2
      NOT2=IROW(LN+1)
*
      IADD10=IAD10(6)

* Long loop, reading buffers until end of buffers is signalled
* by length field holding a negative number.
 300  CONTINUE
      CALL dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
      CALL iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
      LENCOP=ICOP1(nCOP+1)
      IF(LENCOP.EQ.0)GO TO 300
      IF(LENCOP.LT.0)GO TO 350

* Loop over the elements of this buffer
      DO 260 II=1,LENCOP
         INDCOP=ICOP1(II)
         IF(ICHK.NE.0)GO TO 460
         IF(INDCOP.NE.0)GO TO 371
         ICHK=1
         GO TO 260

460      ICHK=0

* Unpack indices NI and NJ from INDCOP
         INDI=INDCOP
*         NI=MOD(INDI,2**10)
*         NJ=MOD(INDI/2**10,2**10)
         NI=IBITS(INDI, 0,10)
         NJ=IBITS(INDI,10,10)

         NSIJ=MUL(NSM(NI),NSM(NJ))
         CALL IPO(IPOF,NVIR,MUL,NSYM,NSIJ,-1)
         IJ1=IROW(NI)+NJ
         ILIM=IPOF(NSYM+1)
* Clear matrices ABIJ, AIBJ, and AJBI.
         CALL VCLR(ABIJ,1,ILIM)
         CALL VCLR(AIBJ,1,ILIM)
         CALL VCLR(AJBI,1,ILIM)
         IF(ITER.EQ.1 .AND. IREST.EQ.0)GO TO 207
*
*     READ (AB/IJ) INTEGRALS
*
         IADR=LASTAD(NOVST+IJ1)
         JTURN=0
 201     CONTINUE

         CALL iDAFILE(Lu_60,2,IBUF,NBITM3+2,IADR)
         CALL dDAFILE(Lu_60,2,BUF,NBITM3,IADR)
         LENBUF=IBUF(NBITM3+1)
         IADR=IBUF(NBITM3+2)
         call faibj5(LENBUF,JTURN,IBUF,BUF, AIBJ,ABIJ)

         IF(IADR.NE.-1) GO TO 201
         IF(JTURN.EQ.1)GO TO 360
*
*     READ (AI/BJ) INTEGRALS
*
207      IADR=LASTAD(NOVST+NOT2+IJ1)
         JTURN=1
         GO TO 201
*
*     CONSTRUCT FIRST ORDER MATRICES
*
360      FAC=1.0D00
         IF(NI.EQ.NJ)FAC=0.5D00
         IN=0
c VV: these calls to getmem are needed to cheat some compilers.

        if (FAC.lt.0) call getmem('CHECK','CHEC','real',0,0)

         IFT=0
         call faibj3(NSIJ,IFT,AIBJ,FSEC,FAC,IN,INS,IPOA,IPOF)


      IF(ITER.EQ.1 .AND. IREST.EQ.0)GO TO 260
      DO 370 IASYM=1,NSYM
        NVIRA=NVIR(IASYM)
        IF(NVIRA.EQ.0)GO TO 370
        IBSYM=MUL(NSIJ,IASYM)
        NVIRB=NVIR(IBSYM)
        IF(NVIRB.EQ.0)GO TO 370
        IPF=IPOF(IASYM)+1
        IPF1=IPOF(IBSYM)+1
        IF(IASYM.GT.IBSYM) THEN
          CALL MTRANS(AIBJ(IPF1),1,AJBI(IPF),1,NVIRA,NVIRB)
          GOTO 370
        END IF
        IF(NSIJ.NE.1) THEN
          CALL MTRANS(ABIJ(IPF1),1,ABIJ(IPF),1,NVIRA,NVIRB)
          CALL MTRANS(AIBJ(IPF1),1,AJBI(IPF),1,NVIRA,NVIRB)
        ELSE
          CALL SQUAR2(ABIJ(IPF),NVIRA)
          IF(NI.EQ.NJ) CALL SQUAR2(AIBJ(IPF),NVIRA)
          CALL MTRANS(AIBJ(IPF),1,AJBI(IPF),1,NVIRA,NVIRB)
        END IF
370   CONTINUE
      GO TO 260
371   CONTINUE
      IF(IFAB.EQ.1) THEN
        CPLA=COP(II)
        IFAB=0
        GO TO 100
      END IF
*      IFAB=MOD(INDCOP,2)
*      ITURN=MOD(INDCOP/2,2)
*      ITYP=MOD(INDCOP/2**2,2**3)
*      ICOUP=MOD(INDCOP/2**5,2**13)
*      ICOUP1=MOD(INDCOP/2**18,2**13)
      IFAB=IBITS(INDCOP, 0, 1)
      ITURN=IBITS(INDCOP, 1, 1)
      ITYP=IBITS(INDCOP, 2, 3)
      ICOUP=IBITS(INDCOP, 5,13)
      ICOUP1=IBITS(INDCOP,18,13)
      CPL=COP(II)
      CPLA=0.0D00
      IF(IFAB.NE.0)GO TO 260
      IF(ITURN.NE.0) GOTO 100
C FIRST ORDER INTERACTION
      INDA=ICOUP
      INDB=IRC(ITYP+1)+ICOUP1
      ISTAR=1
      IF(ITYP.EQ.1)ISTAR=INS+1
      IF(INS.NE.0) THEN
        COPI=CPL*C(INDA)
        CALL DAXPY_(INS,COPI,FSEC(ISTAR),1,S(INDX(INDB)+1),1)
        TERM=DDOT_(INS,FSEC(ISTAR),1,C(INDX(INDB)+1),1)
        S(INDA)=S(INDA)+CPL*TERM
      END IF
      GO TO 260

C INTERACTIONS BETWEEN DOUBLES AND
C INTERACTIONS BETWEEN SINGLES
100   IF((ITER.EQ.1).AND.(IREST.EQ.0))GO TO 260

       Call faibj2(IFTA,IFTB,ICOUP1,ICOUP,
     & INDA,INDB,MYSYM,INTSYM,
     & NYSYM,NSIJ,MYL,NYL,FACS,IPOA,IPOB,
     & INMY,INNY,INDX,iTYP)




      IF(ITYP.NE.5)GO TO 71
C DOUBLET-DOUBLET INTERACTIONS
      IN=IPOF(MYL+1)-IPOF(MYL)
      IF(IN.EQ.0)GO TO 260
      IPF=IPOF(MYL)+1
      CALL DYAX(IN,CPL,AIBJ(IPF),1,F,1)
      CALL DAXPY_(IN,CPLA,ABIJ(IPF),1,F,1)
      IF(INDA.EQ.INDB)CALL SETZZ(F,NVIR(MYL))
*      CALL DGEMTX (NVIR(MYL),NVIR(NYL),FACS,F,NVIR(MYL),
*     *             C(INMY),1,S(INNY),1)
      CALL DGEMV_('T',NVIR(MYL),NVIR(NYL),FACS,F,NVIR(MYL),
     *             C(INMY),1,1.0D0,S(INNY),1)
      IF(INDA.NE.INDB) THEN
*        CALL DGEMX (NVIR(MYL),NVIR(NYL),FACS,F,NVIR(MYL),
*     *              C(INNY),1,S(INMY),1)
        CALL DGEMV_('N',NVIR(MYL),NVIR(NYL),FACS,F,NVIR(MYL),
     *              C(INNY),1,1.0D0,S(INMY),1)
      END IF
      GO TO 260
C TRIPLET-SINGLET, SINGLET-TRIPLET,
C TRIPLET-TRIPLET AND SINGLET-SINGLET INTERACTIONS
71    CONTINUE

      call loop70(INTSYM,INDX,C,S,ABIJ,AIBJ,AJBI,BUF,
     *      IBUF,A,B,F,FSEC,IPOF,IPOA,IPOB,
     * MYL,NYL,INDA,INDB,INMY,INNY,IFTB,IFTA,FACS,
     * IAB,CPL,CPLA, NVIRA,NVIRC,NVIRB)


260   CONTINUE
      GO TO 300

350   CONTINUE
      CALL CSCALE(INDX,INTSYM,C,SQ2INV)
      CALL CSCALE(INDX,INTSYM,S,SQ2)

      RETURN
      END
