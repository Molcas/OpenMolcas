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
      SUBROUTINE INDMAT(ICSPCK,INTSYM,INDX,ISAB,JREFX,CISEL)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "SysDef.fh"

#include "mrci.fh"
      DIMENSION INDX(*),ICSPCK(*),INTSYM(*)
      DIMENSION ISAB(NVIRT,NVIRT),JREFX(*)
      DIMENSION CISEL(NREF,NSEL)
      CHARACTER*20 STR20
CPAM97      EXTERNAL UNPACK
CPAM97      INTEGER UNPACK
CPAM97      JCASE(L)=UNPACK(CSPCK((L+29)/30), 2*L-(2*L-1)/60*60, 2)
      JCASE(L)=ICUNP(ICSPCK,L)
CPAM96      JSYM(L)=UNPACK(INTSYM((L+9)/10),3*MOD(L-1,10)+1,3)+1
      JSYM(L)=JSUNP(INTSYM,L)
C
      ILIM=4
      IF(IFIRST.NE.0)ILIM=2
      DO 5 I=1,NSYM
         NVPAIR(I)=0
5     CONTINUE
      ISMAX=0
      DO 15 NA=1,NVIRT
        DO 25 NB=1,NA
          NSAB=MUL(NSM(LN+NA),NSM(LN+NB))
          NVPAIR(NSAB)=NVPAIR(NSAB)+1
          ISAB(NA,NB)=NVPAIR(NSAB)
          ISMAX=MAX(ISMAX,ISAB(NA,NB))
          ISAB(NB,NA)=ISAB(NA,NB)
25      CONTINUE
        NDIAG(NA)=ISAB(NA,NA)
15    CONTINUE
C INDX - STARTING POINT IN CI VECTOR OF EACH BLOCK WITH A
C COMMON INTERNAL WALK.
C VALENCE CONFIGURATIONS:
      IR=IRC(1)
      DO 10 II=1,IR
        INDX(II)=II
10    CONTINUE
      JSC(1)=IR
C SINGLY EXCITED CONFIGURATIONS:
      IR1=IR+1
      IR2=IRC(2)
      IND=IR
      DO 20 II=IR1,IR2
         INDX(II)=IND
         NSS=MUL(JSYM(II),LSYM)
         IND=IND+NVIR(NSS)
20    CONTINUE
      JSC(2)=IND
      NCDOUB=IND-JSC(1)
      IF(IFIRST.EQ.0) THEN
C DOUBLY EXCITED CONFIGURATIONS:
        IR1=IR2+1
        IR2=IRC(4)
        JSC(3)=JSC(2)
        DO 30 II=IR1,IR2
           INDX(II)=IND
           NSS=MUL(JSYM(II),LSYM)
           IND=IND+NVPAIR(NSS)
           IF(II.EQ.IRC(3))JSC(3)=IND
30      CONTINUE
        JSC(4)=IND
        JJM=(JJS(LSYM+1)-JJS(LSYM))*NVIRT
        NCTRIP=JSC(3)-JSC(2)-JJM
        NCSING=JSC(4)-JSC(3)
      ELSE
        JJM=0
        NCTRIP=0
        NCSING=0
      END IF
      NCONF=JSC(ILIM)
C LIST THE REFERENCE CONFIGURATIONS, AND AT THE SAME TIME,
C IDENTIFY CSFS GIVEN IN SELECTION VECTOR INPUT:
      WRITE(6,*)
      CALL XFLUSH(6)
      WRITE(6,*)'      LIST OF REFERENCE CONFIGURATIONS.'
      CALL XFLUSH(6)
      WRITE(6,*)'     CONF NR:    GUGA CASE NUMBERS OF ACTIVE ORBITALS:'
      CALL XFLUSH(6)
      DO 135 I=1,IRC(1)
        IREF=JREFX(I)
        IF(IREF.EQ.0) GOTO 135
        IREFX(IREF)=I
        IOFF=LN*(I-1)
        WRITE(6,'(5X,I6,7X,30I1)') I,(JCASE(IOFF+J),J=1,LN)
      CALL XFLUSH(6)
        JJ=0
        DO 134 ISEL=1,NSEL
          CISEL(IREF,ISEL)=0.0D00
          NC=NCOMP(ISEL)
          DO 133 IC=1,NC
            STR20=SSEL(JJ+IC)
            DO 132 ILEV=1,LN
              JCAS1=JCASE(IOFF+ILEV)
              READ(STR20(ILEV:ILEV),'(I1)') JCAS2
              IF(JCAS1.NE.JCAS2) GOTO 133
132         CONTINUE
            CISEL(IREF,ISEL)=CSEL(JJ+IC)
            JJ=JJ+NC
            GOTO 134
133       CONTINUE
          JJ=JJ+NC
134     CONTINUE
135   CONTINUE
C ORTHONORMALIZE THE SELECTION VECTORS:
      DO 137 ISEL=1,NSEL
        DO 136 JSEL=1,ISEL-1
          X=DDOT_(NREF,CISEL(1,JSEL),1,CISEL(1,ISEL),1)
          CALL DAXPY_(NREF,-X,CISEL(1,JSEL),1,CISEL(1,ISEL),1)
136     CONTINUE
        X=1.0D00/DDOT_(NREF,CISEL(1,ISEL),1,CISEL(1,ISEL),1)
        CALL DSCAL_(NREF,X,CISEL(1,ISEL),1)
137   CONTINUE
      WRITE(6,*)
      CALL XFLUSH(6)
      WRITE(6,*)'      REAL CONFIGURATIONS:'
      CALL XFLUSH(6)
      IF(IFIRST.EQ.0) THEN
        WRITE(6,215)NREF,NCVAL-NREF,NCDOUB,NCTRIP,NCSING
      CALL XFLUSH(6)
215     FORMAT(/,6X,'               REFERENCE ',I8,
     *         /,6X,'           OTHER VALENCE ',I8,
     *         /,6X,' DOUBLET COUPLED SINGLES ',I8,
     *         /,6X,' TRIPLET COUPLED DOUBLES ',I8,
     *         /,6X,' SINGLET COUPLED DOUBLES ',I8)
      ELSE
        WRITE(6,216)NREF,NCVAL-NREF,NCDOUB
      CALL XFLUSH(6)
216     FORMAT(/,6X,'               REFERENCE ',I8,
     *         /,6X,'           OTHER VALENCE ',I8,
     *         /,6X,' DOUBLET COUPLED SINGLES ',I8)
      END IF
      JSCI=JSC(ILIM)-JJM
      WRITE(6,'(6X,A,I8)')'                  TOTAL ',JSCI
      CALL XFLUSH(6)
      RETURN
      END
