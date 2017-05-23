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
      SUBROUTINE AITD(INTSYM,INDX,C1,C2,TDMO,A,FAK,FKA)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
      DIMENSION INTSYM(*),INDX(*),C1(*),C2(*),
     *          TDMO(NBAST,NBAST),A(*),FAK(*),FKA(*)
      DIMENSION IPOB(9)
*
      JSYM(L)=JSUNP(INTSYM,L)
*
C CALCULATE TRANSITION DENSITY ELEMENTS TDMO(K,A) AND TDMO(A,K),
C WHERE K IS INTERNAL, A IS EXTERNAL ORBITAL.
C SCRATCH AREAS ARE: A(), SIZE NEEDED IS NVMAX**2
C       AND FAK(), FKA(), SIZE NEEDED IS NVMAX
      CALL CSCALE(INDX,INTSYM,C1,SQ2)
      CALL CSCALE(INDX,INTSYM,C2,SQ2)
      NVT=IROW(NVIRT+1)
      ICHK=0
      IJOLD=0
      NK=0
      NSK=1
      IADD10=IAD10(9)
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
      IF(ICHK.NE.0)GO TO 460
      IF(IND.NE.0)GO TO 11
      ICHK=1
      GO TO 10
460   ICHK=0
      IF(IJOLD.NE.0) THEN
C PUT FAK,FKA BACK INTO TDMO.
        NA1=NVIRP(NSK)+1
        NA2=NVIRP(NSK)+NVIR(NSK)
        INK=0
        IF(NA2.LT.NA1)GO TO 10
        DO 113 NA=NA1,NA2
          INK=INK+1
          TDMO(LN+NA,NK)=FAK(INK)
          TDMO(NK,LN+NA)=FKA(INK)
113     CONTINUE
      END IF
      NK=IND
      IJOLD=NK
      NSK=NSM(NK)
C PUT TDMO ELEMENTS INTO ARRAYS FAK, FKA.
      NA1=NVIRP(NSK)+1
      NA2=NVIRP(NSK)+NVIR(NSK)
      INK=0
      IF(NA2.LT.NA1)GO TO 10
      DO 13 NA=NA1,NA2
        INK=INK+1
        FAK(INK)=TDMO(LN+NA,NK)
        FKA(INK)=TDMO(NK,LN+NA)
13    CONTINUE
      GO TO 10
11    IF(INK.EQ.0)GO TO 10
*      ITYP=MOD(IND,2**6)
*      ICP2=MOD(IND/2**6,2**13)
*      ICP1=MOD(IND/2**19,2**13)
      ITYP=IBITS(IND, 0, 6)
      ICP2=IBITS(IND, 6,13)
      ICP1=IBITS(IND,19,13)
      IF(ITYP.GT.1)GO TO 12
      INDA=ICP1
      INDB=IRC(1)+ICP2
      INNY=INDX(INDB)+1
      COPI=C1(INDA)*COP(II)
      CALL DAXPY_(INK,COPI,C2(INNY),1,FAK,1)
      COPI=C2(INDA)*COP(II)
      CALL DAXPY_(INK,COPI,C1(INNY),1,FKA,1)
      GO TO 10
12    CONTINUE
      INDA=IRC(1)+ICP1
      INDB=IRC(ITYP)+ICP2
      INMY=INDX(INDA)+1
      INNY=INDX(INDB)+1
      MYSYM=JSYM(INDA)
      NYSYM=MUL(MYSYM,NSK)
      MYL=MUL(MYSYM,LSYM)
      NYL=MUL(NYSYM,LSYM)
      IFT=0
      IF(ITYP.EQ.2)IFT=1
      CALL IPO(IPOB,NVIR,MUL,NSYM,NYL,IFT)
      NVM=NVIR(MYL)
      COPI=COP(II)
      IF(NYL.NE.1) THEN
        IF(NSK.GT.MYL) THEN
*          CALL DGEMTX (NVM,INK,COPI,C1(INNY+IPOB(NSK)),NVM,
*     *                 C2(INMY),1,FAK,1)
          CALL DGEMV_('T',NVM,INK,COPI,C1(INNY+IPOB(NSK)),NVM,
     *                 C2(INMY),1,1.0D0,FAK,1)
*          CALL DGEMTX (NVM,INK,COPI,C2(INNY+IPOB(NSK)),NVM,
*     *                 C1(INMY),1,FKA,1)
          CALL DGEMV_('T',NVM,INK,COPI,C2(INNY+IPOB(NSK)),NVM,
     *                 C1(INMY),1,1.0D0,FKA,1)
        ELSE
          IF(IFT.EQ.1)COPI=-COPI
*          CALL DGEMX (INK,NVM,COPI,C1(INNY+IPOB(MYL)),INK,
*     *                C2(INMY),1,FAK,1)
          CALL DGEMV_('N',INK,NVM,COPI,C1(INNY+IPOB(MYL)),INK,
     *                C2(INMY),1,1.0D0,FAK,1)
*          CALL DGEMX (INK,NVM,COPI,C2(INNY+IPOB(MYL)),INK,
*     *                C1(INMY),1,FKA,1)
          CALL DGEMV_('N',INK,NVM,COPI,C2(INNY+IPOB(MYL)),INK,
     *                C1(INMY),1,1.0D0,FKA,1)
        END IF
      ELSE
        IF(IFT.EQ.0)CALL SQUAR(C1(INNY+IPOB(MYL)),A,NVM)
        IF(IFT.EQ.1)CALL SQUARN(C1(INNY+IPOB(MYL)),A,NVM)
*        CALL DGEMTX (NVM,INK,COPI,A,NVM,C2(INMY),1,FAK,1)
        CALL DGEMV_('T',NVM,INK,COPI,A,NVM,C2(INMY),1,1.0D0,FAK,1)
        IF(IFT.EQ.0)CALL SQUAR(C2(INNY+IPOB(MYL)),A,NVM)
        IF(IFT.EQ.1)CALL SQUARN(C2(INNY+IPOB(MYL)),A,NVM)
*        CALL DGEMTX (NVM,INK,COPI,A,NVM,C1(INMY),1,FKA,1)
        CALL DGEMV_('T',NVM,INK,COPI,A,NVM,C1(INMY),1,1.0D0,FKA,1)
      END IF
10    CONTINUE
      GO TO 100
200   CONTINUE
      NA1=NVIRP(NSK)+1
      NA2=NVIRP(NSK)+NVIR(NSK)
      INK=0
      DO 213 NA=NA1,NA2
        INK=INK+1
        TDMO(LN+NA,NK)=FAK(INK)
        TDMO(NK,LN+NA)=FKA(INK)
213   CONTINUE
      CALL CSCALE(INDX,INTSYM,C1,SQ2INV)
      CALL CSCALE(INDX,INTSYM,C2,SQ2INV)
      RETURN
      END
