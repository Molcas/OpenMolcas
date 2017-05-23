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
      SUBROUTINE AID(INTSYM,INDX,C,DMO,A,B,FK)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
      DIMENSION INTSYM(*),INDX(*),C(*),DMO(*),
     *          A(*),B(*),FK(*)
      DIMENSION IPOB(9)
*
      JSYM(L)=JSUNP(INTSYM,L)
*
C SCRATCH AREAS: A(),B() AND FK().
      CALL CSCALE(INDX,INTSYM,C,SQ2)
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
C IND=0 INDICATES END OF THIS BLOCK OF COUPLING COEFFS.
      ICHK=1
      GO TO 10
460   CONTINUE
C ICHK=1 INDICATES BEGINNING OF A NEW BLOCK OF COUPLING COEFFS.
      ICHK=0
      IF(IJOLD.NE.0) THEN
C PUT AWAY FK INTO DMO
        NA1=NVIRP(NSK)+1
        NA2=NVIRP(NSK)+NVIR(NSK)
        INK=0
        IF(NA2.LT.NA1)GO TO 10
        DO 113 NA=NA1,NA2
          INK=INK+1
          NAK=IROW(LN+NA)+NK
          DMO(NAK)=FK(INK)
113     CONTINUE
      END IF
      NK=IND
      IJOLD=NK
      NSK=NSM(NK)
C PICK OUT ELEMENTS FROM DMO AND PUT INTO FK:
      NA1=NVIRP(NSK)+1
      NA2=NVIRP(NSK)+NVIR(NSK)
      INK=0
      IF(NA2.LT.NA1)GO TO 10
      DO 13 NA=NA1,NA2
        INK=INK+1
        NAK=IROW(LN+NA)+NK
        FK(INK)=DMO(NAK)
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
      COPI=C(INDA)*COP(II)/ENP
      CALL DAXPY_(INK,COPI,C(INNY),1,FK,1)
      GO TO 10
12    IF(ITER.EQ.1 .AND. IREST.EQ.0)GO TO 10
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
      CALL VCLR(B,1,INK)
      COPI=COP(II)/ENP
      IF(NYL.NE.1) THEN
        IF(NSK.GT.MYL) THEN
          CALL FMMM(C(INMY),C(INNY+IPOB(NSK)),B,1,INK,NVM)
          CALL VSMA(B,1,COPI,FK,1,FK,1,INK)
        ELSE
          CALL FMMM(C(INNY+IPOB(MYL)),C(INMY),B,INK,1,NVM)
          IF(IFT.EQ.1)COPI=-COPI
          CALL VSMA(B,1,COPI,FK,1,FK,1,INK)
        END IF
      ELSE
        IF(IFT.EQ.0)CALL SQUAR(C(INNY+IPOB(MYL)),A,NVM)
        IF(IFT.EQ.1)CALL SQUARN(C(INNY+IPOB(MYL)),A,NVM)
        CALL FMMM(C(INMY),A,B,1,INK,NVM)
        CALL VSMA(B,1,COPI,FK,1,FK,1,INK)
      END IF
10    CONTINUE
      GO TO 100
200   CONTINUE
      NA1=NVIRP(NSK)+1
      NA2=NVIRP(NSK)+NVIR(NSK)
      INK=0
      DO 213 NA=NA1,NA2
        INK=INK+1
        NAK=IROW(LN+NA)+NK
        DMO(NAK)=FK(INK)
213   CONTINUE
      CALL CSCALE(INDX,INTSYM,C,SQ2INV)
      RETURN
      END
