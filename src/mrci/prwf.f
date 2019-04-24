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
      SUBROUTINE PRWF_MRCI(ICSPCK,INTSYM,INDX,C,JREFX)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C(*),INDX(*),ICSPCK(*),
     *          INTSYM(*),JREFX(*)
      CHARACTER*12 CSFTYP
      CHARACTER*14 FORM0,FORM1,FORM2,FORM
      CHARACTER*14 FORM00,FORM01,FORM02

#include "SysDef.fh"

#include "mrci.fh"
      DIMENSION IOC(32),IORBI(32),ISP(32),ILSYM(32)
CPAM97      EXTERNAL UNPACK
CPAM97      INTEGER UNPACK
      DATA FORM00 /'(1X,A,6X,53I2)'/
      DATA FORM01 /'(1X,A,3X,54I2)'/
      DATA FORM02 /'(1X,A,55I2)   '/
      DATA FORM0 /'(1X,A,6X,34I3)'/
      DATA FORM1 /'(1X,A,3X,35I3)'/
      DATA FORM2 /'(1X,A,36I3)   '/
C STATEMENT FUNCTIONS FOR RETRIEVING GUGA CASE NUMBERS AND INTERNAL
C SYMMETRY LABEL:
CPAM97      JO(L)=UNPACK(CSPCK((L+29)/30), 2*L-(2*L-1)/60*60, 2)
      JO(L)=ICUNP(ICSPCK,L)
CPAM96      JSYM(L)=UNPACK(INTSYM((L+9)/10),3*MOD(L-1,10)+1,3)+1
      JSYM(L)=JSUNP(INTSYM,L)
      NA=0
      NB=0
      ILIM=4
      LN2=LN+2
      IF(IFIRST.NE.0)ILIM=2
      SCALE=1.0/SQRT(DDOT_(NCONF,C,1,C,1))
      CALL DSCAL_(NCONF,SCALE,C,1)
      DO 4 J=1,LN
         IORBI(J+2)=IORB(J)
         ILSYM(J+2)=NSM(J)
4     CONTINUE
      JCONF=JSC(1)
      WRITE(6,'(A,F5.3)')'      CI-COEFFICIENTS LARGER THAN ',CTRSH
      CALL XFLUSH(6)
      DO 5 IS=1,NSYM
        IF(NFMO(IS).GT.0) THEN
          WRITE(6,*)' NOTE: THE FOLLOWING ORBITALS WERE FROZEN'
      CALL XFLUSH(6)
          WRITE(6,*)' ALREADY AT THE INTEGRAL TRANSFORMATION STEP'
      CALL XFLUSH(6)
          WRITE(6,*)' AND DO NOT EXPLICITLY APPEAR:'
      CALL XFLUSH(6)
          WRITE(6,'(6X,A,8I4)')'  SYMMETRY:',(I,I=1,NSYM)
      CALL XFLUSH(6)
          WRITE(6,'(6X,A,8I4)')'PRE-FROZEN:',(NFMO(I),I=1,NSYM)
      CALL XFLUSH(6)
          GOTO 6
        END IF
5     CONTINUE
6     CONTINUE
      WRITE(6,*)' ORDER OF SPIN-COUPLING: (PRE-FROZEN, NOT SHOWN)'
      CALL XFLUSH(6)
      WRITE(6,*)'                         (FROZEN, NOT SHOWN)'
      CALL XFLUSH(6)
      WRITE(6,*)'                          VIRTUAL'
      CALL XFLUSH(6)
      WRITE(6,*)'                          ADDED VALENCE'
      CALL XFLUSH(6)
      WRITE(6,*)'                          INACTIVE'
      CALL XFLUSH(6)
      WRITE(6,*)'                          ACTIVE'
      CALL XFLUSH(6)
      WRITE(6,*)
      CALL XFLUSH(6)
      WRITE(6,*)' ORBITALS ARE NUMBERED WITHIN EACH SEPARATE SYMMETRY.'
      CALL XFLUSH(6)
      DO 10 I=1,NCONF
        CI=C(I)
        ACI=ABS(C(I))
        IF(I.LE.JCONF) THEN
          JMIN=1
          NREXT=0
          IF(JREFX(I).NE.0) THEN
            CSFTYP='   REFERENCE'
            CLIM=0.0D00
          ELSE
            CSFTYP='     VALENCE'
            CLIM=CTRSH
          END IF
        ELSE IF (I.LE.JSC(2)) THEN
          NREXT=1
          JMIN=1+IRC(1)
          CSFTYP='     DOUBLET'
          CLIM=CTRSH
        ELSE IF (I.LE.JSC(3)) THEN
          NREXT=2
          JMIN=1+IRC(2)
          CSFTYP='     TRIPLET'
          CLIM=CTRSH
        ELSE
          NREXT=2
          JMIN=1+IRC(3)
          CSFTYP='     SINGLET'
          CLIM=CTRSH
        END IF
        IF(ACI.LT.CLIM) GOTO 10
         JJ=I
         IJ=I
         IF(NREXT.GT.0) THEN
           JMAX=IRC(ILIM)
           DO 20 J=JMIN,JMAX
              JJ=J
              IF(INDX(J).LT.IJ)GO TO 20
              JJ=JJ-1
              GOTO 25
20         CONTINUE
25         CONTINUE
         END IF
         NSJ=MUL(JSYM(JJ),LSYM)
         JVIR=IJ-INDX(JJ)
         II1=(JJ-1)*LN
         DO 31 II=1,LN
            II1=II1+1
            ISP(II+2)=JO(II1)
            IOC(II+2)=(1+ISP(II+2))/2
31       CONTINUE
         IF(NREXT.EQ.0) THEN
           FORM=FORM0
           IF(LN2.GT.36) FORM=FORM00
           LN1=3
         ELSE IF (NREXT.EQ.1) THEN
           IO=JVIR+NVIRP(NSJ)+LN
           IORBI(2)=IORB(IO)
           IOC(2)=1
           ISP(2)=1
           ILSYM(2)=NSJ
           FORM=FORM1
           IF(LN2.GT.36) FORM=FORM01
           LN1=2
         ELSE
           IN=0
           DO 46 II=1,NVIRT
              NA=II
              NSI=MUL(NSJ,NSM(LN+II))
              J1=NVIRP(NSI)+1
              J2=NVIRP(NSI)+NVIR(NSI)
              IF(J2.GT.II)J2=II
              IF(J2.LT.J1)GO TO 46
              DO 47 J=J1,J2
                 NB=J
                 IN=IN+1
                 IF(IN.EQ.JVIR)GO TO 48
47            CONTINUE
46         CONTINUE
48         CONTINUE
           IOC(1)=1
           ISP(1)=1
           ILSYM(1)=NSM(LN+NB)
           IO=LN+NB
           IORBI(1)=IORB(IO)
           IF(NA.EQ.NB) THEN
             IORBI(2)=IORBI(1)
             IOC(2)=2
             ISP(2)=3
             ILSYM(2)=NSM(IO)
             CI=CI*SQRT(0.5D00)
             IF(ABS(CI).LT.CTRSH)GO TO 10
             FORM=FORM1
             IF(LN2.GT.36) FORM=FORM01
             LN1=2
           ELSE
             IOC(2)=1
             ISP(2)=2
             IF(CSFTYP.EQ.'     TRIPLET') ISP(2)=1
             IO=LN+NA
             IORBI(2)=IORB(IO)
             ILSYM(2)=NSM(IO)
             FORM=FORM2
             IF(LN2.GT.36) FORM=FORM02
             LN1=1
           END IF
         END IF
         WRITE(6,*)
      CALL XFLUSH(6)
         WRITE(6,105) I,C(I),CSFTYP
      CALL XFLUSH(6)
105   FORMAT(/6X,'CONFIGURATION',I7,3X,'COEFFICIENT',F10.6,A)
         WRITE(6,FORM) 'SYMMETRY     ',(ILSYM(J),J=LN1,LN2)
      CALL XFLUSH(6)
         WRITE(6,FORM) 'ORBITALS     ',(IORBI(J), J=LN1,LN2)
      CALL XFLUSH(6)
         WRITE(6,FORM) 'OCCUPATION   ',(IOC(J),  J=LN1,LN2)
      CALL XFLUSH(6)
         WRITE(6,FORM) 'SPIN-COUPLING',(ISP(J),  J=LN1,LN2)
      CALL XFLUSH(6)
10    CONTINUE
      RETURN
      END
