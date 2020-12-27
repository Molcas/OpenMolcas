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
*      SUBROUTINE AI(INTSYM,INDX,C,S,FC,BUFIN,IBUFIN,A,B,FK,DBK,KTYP)
*      SUBROUTINE AI_MRCI(INTSYM,INDX,C,S,FC,BUF,IBUF,A,B,FK,DBK,KTYP)
      SUBROUTINE AI_MRCI(INTSYM,INDX,C,S,FC,A,B,FK,DBK,KTYP)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "mrci.fh"
      DIMENSION INTSYM(*),INDX(*),C(*),S(*),
*     *          FC(*),BUFIN(*),IBUFIN(*),
*     *          FC(*),BUF(NBITM3),IBUF(NBITM3+2),
     *          FC(*),
     *          A(*),B(*),FK(*),DBK(*)
      DIMENSION IPOB(9)
      PARAMETER (ONE=1.0D00)
*
      JSYM(L)=JSUNP(INTSYM,L)
*
C KTYP=0,  (A/I)   INTEGRALS
C KTYP=1,  (AI/JK) INTEGRALS

      CALL GETMEM('BUF','ALLO','REAL',LBUF,NBITM3)
      CALL GETMEM('IBUF','ALLO','INTE',LIBUF,NBITM3+2)

      CALL CSCALE(INDX,INTSYM,C,SQ2)
      CALL CSCALE(INDX,INTSYM,S,SQ2INV)
      NVT=IROW(NVIRT+1)
      ICHK=0
      IJOLD=0
      NK=0
      NSA=1
      NOTT=LN*(LN+1)
      NOVST=LN*NVIRT+1+NVT
CPAM97 New portable code:
      NBCMX3=(RTOI*NBSIZ3-2)/(RTOI+1)
      IBOFF3=RTOI*NBCMX3
      IBBC3=IBOFF3+NBCMX3+1
*PAM04      IBDA3=IBBC3+1

      IF(KTYP.EQ.0)IADD10=IAD10(9)
      IF(KTYP.EQ.1)IADD10=IAD10(7)
C READ A COP BUFFER
100   CONTINUE
      CALL dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
      CALL iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
      LEN=ICOP1(nCOP+1)
      IF(LEN.EQ.0)GO TO 100
      IF(LEN.LT.0)GO TO 200
C LOOP THROUGH THE COP BUFFER:
      DO 10 II=1,LEN
      IND=ICOP1(II)
      IF(ICHK.NE.0) THEN
C BEGIN A RATHER LONG IF-BLOCK.
C ICHK FLAG IS SET. THIS SIGNALS THAT PREVIOUS IND WAS 0, WHICH IS
C USED TO INDICATE CHANGE TO A NEW BLOCK OF COUPLING COEFFICIENTS.
C RESET ICHK FLAG.
        ICHK=0
        IF(KTYP.EQ.0) THEN
C AI CASE. SAVE INTERNAL ORBITAL INDEX IN NK:
          NK=IND
          IJOLD=NK
          NSK=NSM(NK)
          NSA=NSK
          GO TO 20
        END IF
C AIJK CASE. UNPACK INTERNAL ORBITAL INDICES INTO NI,NJ,NK:
        INDI=IND
*        NI=MOD(INDI,2**10)
*        NJ=MOD(INDI/2**10,2**10)
*        NK=MOD(INDI/2**20,2**10)
        NI=IBITS(INDI, 0,10)
        NJ=IBITS(INDI,10,10)
        NK=IBITS(INDI,20,10)
        NSI=NSM(NI)
        NSJ=NSM(NJ)
        NSK=NSM(NK)
        NSIJ=MUL(NSI,NSJ)
        NSA=MUL(NSIJ,NSK)
        IJ=IROW(NI)+NJ
        IF(IJ.NE.IJOLD) THEN
C NEW INTERNAL PAIR IJ. LOAD A NEW SET OF INTEGRALS INTO FC:
          IJOLD=IJ
          IADR=LASTAD(NOVST+NOTT+IJ)
          CALL FZERO(FC,NBTRI)

90        CONTINUE
*PAM04          CALL dDAFILE(Lu_60,2,IBUFIN,NBSIZ3,IADR)
          CALL iDAFILE(Lu_60,2,iWORK(LIBUF),NBITM3+2,IADR)
          CALL dDAFILE(Lu_60,2,WORK(LBUF),NBITM3,IADR)
          LENGTH=iWORK(LIBUF+NBITM3)
          IADR  =iWORK(LIBUF+NBITM3+1)
*PAM04          LENGTH=IBUFIN(IBBC3)
*PAM04          IADR=IBUFIN(IBDA3)
          IF(LENGTH.EQ.0)GO TO 91
*          CALL SCATTER(LENGTH,FC,IBUFIN(IBOFF3+1),BUFIN)
          do i=0,length-1
*PAM04            fc(IBUFIN(IBOFF3+i))=bufin(i)
            fc(iWORK(LIBUF+i))=WORK(LBUF+i)
          end do
91        IF(IADR.NE.-1) GO TO 90
        END IF
20      CONTINUE
C FOR THIS PARTICULAR K, TRANSFER FC(NK,NA) TO ARRAY FK:
        NVIRA=NVIR(NSA)
        IF(NVIRA.EQ.0) GOTO 10
        DO 13 I=1,NVIRA
          NA=NVIRP(NSA)+I
          NAK=IROW(LN+NA)+NK
          FK(I)=FC(NAK)
13      CONTINUE
        GOTO 10
      END IF
C END OF THE LONG IF-BLOCK.
      IF(IND.EQ.0) THEN
C IND=0 SIGNALS SWITCH TO A NEW SET OF INTEGRALS.
        ICHK=1
        GO TO 10
      END IF
C WE ARE PROCESSING A COUPLING COEFFICIENT AS USUAL.
      IF(NVIRA.EQ.0)GO TO 10
*      ITYP=MOD(IND,2**6)
*      ICP2=MOD(IND/2**6,2**13)
*      ICP1=MOD(IND/2**19,2**13)
      ITYP=IBITS(IND, 0, 6)
      ICP2=IBITS(IND, 6,13)
      ICP1=IBITS(IND,19,13)
      IF(ITYP.GT.1)GO TO 12
C ITYP=1. VALENCE-SINGLES CASE.
      INDA=ICP1
      INDB=IRC(1)+ICP2
      INNY=INDX(INDB)+1
      COPI=COP(II)*C(INDA)
      CALL DAXPY_(NVIRA,COPI,FK,1,S(INNY),1)
      TERM=DDOT_(NVIRA,FK,1,C(INNY),1)
      S(INDA)=S(INDA)+COP(II)*TERM
      GO TO 10
12    IF(ITER.EQ.1 .AND. IREST.EQ.0)GO TO 10
      INDA=IRC(1)+ICP1
      INDB=IRC(ITYP)+ICP2
      INMY=INDX(INDA)+1
      INNY=INDX(INDB)+1
      MYINTS=JSYM(INDA)
      NYINTS=MUL(MYINTS,NSA)
      MYEXTS=MUL(MYINTS,LSYM)
      NYEXTS=MUL(NYINTS,LSYM)
      IFT=0
      IF(ITYP.EQ.2)IFT=1
      CALL IPO(IPOB,NVIR,MUL,NSYM,NYEXTS,IFT)
      NVM=NVIR(MYEXTS)
      CALL FZERO(DBK,NVIRA)
      CALL DAXPY_(NVIRA,COP(II),FK,1,DBK,1)
      IF(NYEXTS.NE.1)GO TO 25
      IF(IFT.EQ.0)CALL SQUAR(C(INNY+IPOB(MYEXTS)),A,NVM)
      IF(IFT.EQ.1)CALL SQUARM(C(INNY+IPOB(MYEXTS)),A,NVM)
      CALL FZERO(B,NVM)
      CALL FMMM(DBK,A,B,1,NVM,NVIRA)
      CALL DAXPY_(NVM,ONE,B,1,S(INMY),1)
      SIGN=1.0D00
      IF(IFT.EQ.1)SIGN=-1.0D00
      IOUT=INNY+IPOB(MYEXTS)-1
      DO 125 I=1,NVM
      DO 130 J=1,I
      IOUT=IOUT+1
      TERM=DBK(I)*C(INMY+J-1)+SIGN*DBK(J)*C(INMY+I-1)
      S(IOUT)=S(IOUT)+TERM
130   CONTINUE
      IF(IFT.EQ.1)GO TO 125
      TERM=DBK(I)*C(INMY+I-1)
      S(IOUT)=S(IOUT)-TERM
125   CONTINUE
      GO TO 10
25    NKM=NVIRA*NVM
      CALL FZERO(B,NVM)
      IF(NSA.GT.MYEXTS)GO TO 26
      IF(IFT.EQ.1)CALL VNEG(DBK,1,DBK,1,NVIRA)
      CALL FMMM(DBK,C(INNY+IPOB(MYEXTS)),B,1,NVM,NVIRA)
      CALL DAXPY_(NVM,ONE,B,1,S(INMY),1)
      CALL FZERO(B,NKM)
      CALL FMMM(DBK,C(INMY),B,NVIRA,NVM,1)
      CALL DAXPY_(NKM,ONE,B,1,S(INNY+IPOB(MYEXTS)),1)
      GO TO 10
26    CALL FMMM(C(INNY+IPOB(NSA)),DBK,B,NVM,1,NVIRA)
      CALL DAXPY_(NVM,ONE,B,1,S(INMY),1)
      CALL FZERO(B,NKM)
      CALL FMMM(C(INMY),DBK,B,NVM,NVIRA,1)
      CALL DAXPY_(NKM,ONE,B,1,S(INNY+IPOB(NSA)),1)
      GO TO 10
10    CONTINUE
      GO TO 100
200   CONTINUE
      CALL CSCALE(INDX,INTSYM,C,SQ2INV)
      CALL CSCALE(INDX,INTSYM,S,SQ2)
      CALL GETMEM('BUF','FREE','REAL',LBUF,NBITM3)
      CALL GETMEM('IBUF','FREE','INTE',LIBUF,NBITM3+2)
      RETURN
      END
