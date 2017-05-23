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
      SUBROUTINE NATORB(D,CM,CMO,DSYM,CAO,OCC,M)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION D(*),CM(*),CMO(*),DSYM(*),CAO(*),OCC(*)

#include "SysDef.fh"

#include "cpfmcpf.fh"
C     FORM DENSITY MATRIX BY SYMMETRY
      NBM=NBAS(M)
      IF(NBM.EQ.0)RETURN
      IF(NORB(M).EQ.0)RETURN
      IF(IPRINT.GE.15)WRITE(6,1)M
1     FORMAT(///,5X,'SYMMETRY NUMBER',I3)
      NBP=0
      M1=M-1
      DO 7 I=1,M1
      NBP=NBP+NORB(I)
7     CONTINUE
      NFR=NFRO(M)
      NORBM=NFR+NISH(M)+NASH(M)+NVIR(M)
      IF(NORBM.EQ.0) Return
      NORBM2=IROW(NORBM+1)
      DO 5 I=1,NORBM2
         DSYM(I)=0.0D0
5     CONTINUE
      IF(NFR.EQ.0)GO TO 10
      II=0
      DO 11 I=1,NFR
         II=II+I
         DSYM(II)=2.0D0
11    CONTINUE
C     REST OF DENSITY MATRIX
10    IJ=0
      DO 15 I=1,NORBM
         DO 20 J=1,I
            IJ=IJ+1
            NI=ICH(NBP+I)
            IF(NI.LT.0)GO TO 20
            NJ=ICH(NBP+J)
            IF(NJ.LT.0)GO TO 20
            NIJ=IROW(NI)+NJ
            IF(NJ.GT.NI)NIJ=IROW(NJ)+NI
            DSYM(IJ)=D(NIJ)
20       CONTINUE
15    CONTINUE
C     DIAGONALIZE
      CALL JACSCF(DSYM,CMO,OCC,NORBM,-1,1.D-11)
      DO 80 I=1,NORBM
      OCC(I)=-OCC(I)
80    CONTINUE
      CALL ORDER(CMO,OCC,NORBM)
      IF(IPRINT.GE.15)WRITE(6,30)
30    FORMAT(//,5X,'NATURAL ORBITALS IN MO-BASIS',//,
     *7X,'OCCUPATION NUMBER',5X,'COEFFICIENTS')
      ILAS=0
      DO 35 I=1,NORBM
      ISTA=ILAS+1
      ILAS=ILAS+NORBM
      OCC(I)=-OCC(I)
      IF(IPRINT.GE.15)WRITE(6,40)I,OCC(I),(CMO(J),J=ISTA,ILAS)
40    FORMAT(/,5X,I4,F10.6,5F10.6,/(19X,5F10.6))
35    CONTINUE
C     TRANSFORM TO AO-BASIS
      IF(IPRINT.GE.15)WRITE(6,45)
45    FORMAT(//,5X,'NATURAL ORBITALS IN AO-BASIS',//,
     *11X,'OCCUPATION NUMBER',5X,'COEFFICIENTS')
      IJ0=-NORBM
      kp = 1
      DO 50 I=1,NORBM
         IJ0=IJ0+NORBM
*
         DO 60 IP=1,NBM
            TERM=D0
            IJ=IJ0
            JP=IP-NBM
            DO 70 J=1,NORBM
               IJ=IJ+1
               JP=JP+NBM
               TERM=TERM+CMO(IJ)*CM(JP)
70          CONTINUE
            CAO(kp+IP-1)=TERM
60       CONTINUE
*
         IF(IPRINT.GE.15)WRITE(6,40)I,OCC(I),(CAO(IP),IP=kp,kp+nbm-1)
         kp = kp + NBM
50    CONTINUE
      RETURN
      END
