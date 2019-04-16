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
      SUBROUTINE INDMAT_CPF(JSY,INDEX,ISAB,ISMAX,JREFX)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "SysDef.fh"

#include "cpfmcpf.fh"
      DIMENSION INDEX(*),JSY(*),ISAB(*),JREFX(*)
      DIMENSION ICOUS(8)
CPAM97      EXTERNAL UNPACK
CPAM97      INTEGER UNPACK
CRL   JSYM(L)=IAND(ISHFT(JSY((L+19)/20),-3*((L+19)/20*20-L)),7)+1
CPAM96      JSYM(L)=UNPACK(JSY((L+9)/10),3*MOD(L-1,10)+1,3)+1
      JSYM(L)=JSUNP_CPF(JSY,L)
C
C     DETERMINE REFERENCE STATE
      JCONF=ISC(1)
      DO 31 IR=1,JCONF
      IF(JREFX(IR).EQ.1)IREF0=IR
31    CONTINUE
      IF(IPRINT.GT.5)WRITE(6,999)IREF0,(JREFX(IR),IR=1,JCONF)
999   FORMAT(2X,I3,2X,'JREFX',10I5)
C
      ILIM=4
      IF(IFIRST.NE.0)ILIM=2
      NSYS(1)=0
      IF(NSYM.EQ.1)GO TO 1
      DO 2 I=2,NSYM
      NSYS(I)=NSYS(I-1)+NVIR(I-1)
2     CONTINUE
1     NSYS(NSYM+1)=NVIRT
      DO 5 I=1,NSYM
      ICOUS(I)=0
      NNS(I)=0
5     CONTINUE
      ISMAX=0
      IN0=-NVIRT
      DO 15 NA=1,NVIRT
      IN0=IN0+NVIRT
      IN=IN0
      IN2=-NVIRT+NA
      DO 25 NB=1,NA
      IN=IN+1
      IN2=IN2+NVIRT
      NSAB=MUL(NSM(LN+NA),NSM(LN+NB))
      ICOUS(NSAB)=ICOUS(NSAB)+1
      ISAB(IN)=ICOUS(NSAB)
      IF(ISMAX.LT.ISAB(IN))ISMAX=ISAB(IN)
      ISAB(IN2)=ISAB(IN)
      IF(ISAB(IN).GT.NNS(NSAB))NNS(NSAB)=ISAB(IN)
25    CONTINUE
      NDIAG(NA)=ISAB(IN)
15    CONTINUE
      IND=0
      IR=IRC(1)
      DO 10 II=1,IR
      IND=IND+1
      INDEX(II)=IND
10    CONTINUE
      JSC(1)=IND
      IR1=IR+1
      IR2=IRC(2)
      DO 20 II=IR1,IR2
      INDEX(II)=IND
      NSS=MUL(JSYM(II),LSYM)
      IND=IND+NVIR(NSS)
20    CONTINUE
      JSC(2)=IND
      IF(IFIRST.NE.0)GO TO 22
      IR1=IR2+1
      IR2=IRC(4)
      JSC(3)=IND
      DO 30 II=IR1,IR2
      INDEX(II)=IND
      NSS=MUL(JSYM(II),LSYM)
      IND=IND+ICOUS(NSS)
      IF(II.EQ.IRC(3))JSC(3)=IND
30    CONTINUE
      JSC(4)=IND
22    IX1=JSC(1)
      IX2=JSC(2)-JSC(1)
      WRITE(6,213)
      CALL XFLUSH(6)
213   FORMAT(//,6X,'FULL-SPACE CONFIGURATIONS (REAL)')
      IF(IFIRST.NE.0)GO TO 212
      JJM=(JJS(LSYM+1)-JJS(LSYM))*NVIRT
      IX3=JSC(3)-JSC(2)-JJM
      IX4=JSC(4)-JSC(3)
      WRITE(6,215)IX1,IX2,IX3,IX4
      CALL XFLUSH(6)
215   FORMAT(/,6X,'NUMBER OF VALENCE STATES',I16,
     */,6X,'NUMBER OF DOUBLET COUPLED SINGLES',I7,
     */,6X,'NUMBER OF TRIPLET COUPLED DOUBLES',I7,
     */,6X,'NUMBER OF SINGLET COUPLED DOUBLES',I7)
      GO TO 211
212   WRITE(6,216)IX1,IX2
      CALL XFLUSH(6)
216   FORMAT(/,6X,'NUMBER OF VALENCE STATES',I14,
     */,6X,'NUMBER OF DOUBLET COUPLED SINGLES',I7)
      JJM=0
211   JSCI=JSC(ILIM)-JJM
      WRITE(6,50)ISC(ILIM),JSCI
      CALL XFLUSH(6)
50    FORMAT(//6X,'FORMAL NUMBER OF CONFIGURATIONS',I8,
     */8X,'REAL NUMBER OF CONFIGURATIONS',I8)
      RETURN
      END
