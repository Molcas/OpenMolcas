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
      SUBROUTINE PRWF(ICASE,JSY,INDEX,C)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C(*),INDEX(*),JSY(*)
      DIMENSION ICASE(*)

#include "SysDef.fh"

#include "cpfmcpf.fh"
      DIMENSION IOC(57),IORB(57),ISP(57),ILSYM(57)
CPAM97      EXTERNAL UNPACK
CPAM97      INTEGER UNPACK
CPAM97      JO(L)=UNPACK(QOCC((L+29)/30), 2*L-(2*L-1)/60*60, 2)
      JO(L)=ICUNP(ICASE,L)
CPAM96      JSYM(L)=UNPACK(JSY((L+9)/10),3*MOD(L-1,10)+1,3)+1
      JSYM(L)=JSUNP(JSY,L)
      NA = 0 ! dummy initialized
      NB = 0 ! dummy initialized
      ILIM=4
      IF(IFIRST.NE.0)ILIM=2
      NCONF=JSC(ILIM)
      IF(ISDCI.EQ.1) CALL DSCAL_(NCONF,1.0D0/DNORM2(NCONF,C,1),C,1)
      JCONF=JSC(1)
      THRC=CTRSH
      WRITE(6,5)THRC
      CALL XFLUSH(6)
      IF(ISDCI.EQ.0)WRITE(6,6)
*
      DO 4 J=1,LN
         IORB(J+2)=J
         ILSYM(J+2)=NSM(J)
4     CONTINUE
*
      DO 10 I=1,NCONF
         JJ=I
         IJ=I
         IF (I.EQ.IREF0) THEN
            WRITE(6,105)I,C(I),'REFERENCE'
      CALL XFLUSH(6)
            GO TO 26
         END IF
         CI=C(I)
         IF (ABS(CI).LT.THRC)GO TO 10
         IF (I.LE.JCONF) THEN
            WRITE(6,105)I,CI,'VALANCE'
      CALL XFLUSH(6)
            GO TO 26
         END IF
         IF (I.LE.JSC(2)) THEN
            JMIN=IRC(1)+1
            WRITE(6,105)I,CI,'DOUBLET'
      CALL XFLUSH(6)
         ELSE IF (I.LE.JSC(3)) THEN
            JMIN=IRC(2)+1
         ELSE
            JMIN=IRC(3)+1
         END IF
         IX1=IRC(ILIM)
         DO 20 J=JMIN,IX1
            JJ=J-1
            IF (INDEX(J).GE.IJ)GO TO 25
20       CONTINUE
25       CONTINUE
26       CONTINUE
         NSJ=MUL(JSYM(JJ),LSYM)
         JVIR=I-INDEX(JJ)
         IF (I.GT.JCONF)JVIR=IJ-INDEX(JJ)
         II1=(JJ-1)*LN
         DO 31 II=1,LN
            II1=II1+1
            ISP(II+2)=JO(II1)
            JOJ=ISP(II+2)
            IF (JOJ.GT.1)JOJ=JOJ-1
            IOC(II+2)=JOJ
31       CONTINUE
         IF (JJ.LE.IRC(1)) THEN
            IORB(1)=0
            IOC(1)=0
            ISP(1)=0
            ILSYM(1)=0
            IORB(2)=0
            IOC(2)=0
            ISP(2)=0
            ILSYM(2)=0
            GO TO 100
         END IF
         IF (JJ.LE.IRC(2)) THEN
            IORB(2)=JVIR+NSYS(NSJ)+LN
            IOC(2)=1
            ISP(2)=1
            ILSYM(2)=NSJ
            IORB(1)=0
            IOC(1)=0
            ISP(1)=0
            ILSYM(1)=0
            GO TO 100
         END IF
         IN=0
         DO 46 II=1,NVIRT
            NA=II
            NSI=MUL(NSJ,NSM(LN+II))
            J1=NSYS(NSI)+1
            J2=NSYS(NSI+1)
            IF (J2.GT.II)J2=II
            IF (J2.LT.J1)GO TO 46
            DO 47 J=J1,J2
               NB=J
               IN=IN+1
               IF (IN.EQ.JVIR)GO TO 48
47          CONTINUE
46       CONTINUE
48       CONTINUE
         IORB(1)=LN+NB
         IOC(1)=1
         ISP(1)=1
         ILSYM(1)=NSM(IORB(1))
         IF (NA.EQ.NB) THEN
            IORB(2)=IORB(1)
            IOC(2)=2
            ISP(2)=3
            ILSYM(2)=NSM(IORB(2))
            IORB(1)=0
            IOC(1)=0
            ISP(1)=0
            ILSYM(1)=0
            CI=CI/SQ2
            IF (ABS(CI).LT.THRC)GO TO 10
         ELSE
            IORB(2)=LN+NA
            IOC(2)=1
            ISP(2)=2
            ILSYM(2)=NSM(IORB(2))
         END IF
         IF (JJ.LE.IRC(3))WRITE(6,105)I,CI,'TRIPLET'
         IF (JJ.GT.IRC(3))WRITE(6,105)I,CI,'SINGLET'
100      CONTINUE
         WRITE(6,*)
      CALL XFLUSH(6)
         IF(LN+2.LE.36)THEN
          WRITE(6,120) 'ORBITALS     ',(IORB(J), J=1,LN+2)
      CALL XFLUSH(6)
          WRITE(6,120) 'OCCUPATION   ',(IOC(J),  J=1,LN+2)
      CALL XFLUSH(6)
          WRITE(6,120) 'SPIN-COUPLING',(ISP(J),  J=1,LN+2)
      CALL XFLUSH(6)
          WRITE(6,120) 'SYMMETRY     ',(ILSYM(J),J=1,LN+2)
      CALL XFLUSH(6)
         ELSE
          WRITE(6,121) 'ORBITALS     ',(IORB(J), J=1,LN+2)
      CALL XFLUSH(6)
          WRITE(6,121) 'OCCUPATION   ',(IOC(J),  J=1,LN+2)
      CALL XFLUSH(6)
          WRITE(6,121) 'SPIN-COUPLING',(ISP(J),  J=1,LN+2)
      CALL XFLUSH(6)
          WRITE(6,121) 'SYMMETRY     ',(ILSYM(J),J=1,LN+2)
      CALL XFLUSH(6)
         END IF
10    CONTINUE
      RETURN
5     FORMAT(//6X,'PRINTOUT OF CI-COEFFICIENTS LARGER THAN',F10.2,/)
6     FORMAT(/6X,'WAVE FUNCTION NOT NORMALIZED',/)
105   FORMAT(/6X,'CONFIGURATION',I7,3X,'COEFFICIENT',F10.6,3X,A)
120   FORMAT(6X,A,36I3)
121   FORMAT(6X,A,55I2)
      END
