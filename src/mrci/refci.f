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
      SUBROUTINE REFCI(HREF,AREF,EREF,ICSPCK,CISEL,PLEN)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "SysDef.fh"

#include "mrci.fh"
      DIMENSION HREF((NREF*(NREF+1))/2),AREF(NREF,NREF),EREF(NREF)
      DIMENSION PLEN(NREF),CISEL(NREF,NSEL),ICSPCK(NCSPCK)
      CHARACTER*2 STR
      CHARACTER*48 FORM1,FORM2,FORM3,FORM4
CPAM97      EXTERNAL UNPACK
CPAM97      INTEGER UNPACK
CPAM97      JCASE(L)=UNPACK(CSPCK((L+29)/30), 2*L-(2*L-1)/60*60, 2)
      JCASE(L)=ICUNP(ICSPCK,L)
      CALL QENTER('REFCI')
      WRITE(6,*)
      CALL XFLUSH(6)
      WRITE(6,*)('-',I=1,60)
      CALL XFLUSH(6)
      WRITE(6,*)'   REFERENCE CI CALCULATION.'
      CALL XFLUSH(6)
      WRITE(6,*)('-',I=1,60)
      CALL XFLUSH(6)
      IF(NSEL.EQ.0) THEN
        WRITE(6,*)' ROOT SELECTION BY ENERGY ORDERING.'
      CALL XFLUSH(6)
        IF(NRROOT.EQ.1) THEN
          WRITE(6,'(A,I8)')'  ONE SINGLE ROOT, NUMBER.....: ',IROOT(1)
      CALL XFLUSH(6)
        ELSE
          WRITE(6,*)' THE FOLLOWING ROOTS WILL BE SELECTED:'
      CALL XFLUSH(6)
          WRITE(6,'(12(A,I2))') ' ROOTS NR ',IROOT(1),
     *                           (',',IROOT(I),I=2,NRROOT-1),
     *                          ', AND ',IROOT(NRROOT)
      CALL XFLUSH(6)
        END IF
      ELSE
        WRITE(6,*)' ROOT SELECTION BY PROJECTION: THE EIGENVECTORS OF'
      CALL XFLUSH(6)
        WRITE(6,*)' THE REFERENCE CI ARE ORDERED BY DECREASING SIZE OF'
      CALL XFLUSH(6)
        WRITE(6,*)' THEIR PROJECTIONS ONTO A SELECTION SPACE.'
      CALL XFLUSH(6)
        IF(NRROOT.EQ.1) THEN
          WRITE(6,*)' SELECT THE EIGENVECTOR WITH LARGEST PROJECTION.'
      CALL XFLUSH(6)
        ELSE
          WRITE(6,'(A,I2,A)')' SELECT THE ',NRROOT,
     *              ' EIGENVECTORS WITH LARGEST PROJECTION.'
      CALL XFLUSH(6)
        END IF
        WRITE(6,*)' THE SELECTION SPACE IS SPANNED BY THE FOLLOWING',
     *            ' VECTORS (NONZERO COMPONENTS ONLY):'
      CALL XFLUSH(6)
        JJ=0
        DO 1234 I=1,NSEL
          WRITE(6,'(A,I2)') ' VECTOR NR. ',I
      CALL XFLUSH(6)
          NC=NCOMP(I)
          WRITE(6,'(5X,I2,5X,A20,F12.8)')
     *          (J,SSEL(JJ+J),CSEL(JJ+J),J=1,NC)
      CALL XFLUSH(6)
          JJ=JJ+NC
1234  CONTINUE
      END IF
      WRITE(6,*)
      CALL XFLUSH(6)
      CALL JACSCF(HREF,AREF,EREF,NREF,-1,1.0D-11)
      CALL ORDER(AREF,EREF,NREF)
      CALL CI_SELECT_MRCI(NREF,AREF,PLEN,NSEL,CISEL,NRROOT,IROOT)
      IF(NSEL.GT.0) THEN
          WRITE(6,*)' THE FOLLOWING ROOTS WERE SELECTED:'
      CALL XFLUSH(6)
          WRITE(6,'(12(A,I2))') ' ROOTS NR ',IROOT(1),
     *                           (',',IROOT(I),I=2,NRROOT-1),
     *                          ', AND ',IROOT(NRROOT)
      CALL XFLUSH(6)
      END IF
      WRITE(STR,'(I2)') LN
      CALL XFLUSH(6)
      FORM1='(2X,'//STR//'X,A,I7,2(8X,I7))'
      FORM2='(2X,'//STR//'X,A,3F15.8)'
      FORM3='('' CSF NR'',I5,'' CASE '','//STR//'I1,3(F13.6,2X))'
      FORM4='(''       '',I5,''      '','//STR//'I1,3(F13.6,2X))'
      WRITE(6,*)
      CALL XFLUSH(6)
      WRITE(6,*)'        LOWEST REFERENCE CI ROOTS:'
      CALL XFLUSH(6)
      NPRT=MIN(NREF,IROOT(NRROOT)+2)
      DO 100 K1=1,NPRT,3
        K2=MIN(NPRT,K1+2)
        WRITE(6,FORM1)'            ROOT',(K,K=K1,K2)
      CALL XFLUSH(6)
        IF(NSEL.GT.0) WRITE(6,FORM2)'SELECTION WEIGHT',(PLEN(K),K=K1,K2)
        WRITE(6,FORM2)'          ENERGY',(EREF(K),K=K1,K2)
      CALL XFLUSH(6)
        DO  50 IREF=1,NREF
          IC=IREFX(IREF)
          IOFF=LN*(IC-1)
          IF(IREF.EQ.1) THEN
          WRITE(6,FORM3)IC,(JCASE(IOFF+J),J=1,LN),(AREF(IREF,K),K=K1,K2)
      CALL XFLUSH(6)
          ELSE
          WRITE(6,FORM4)IC,(JCASE(IOFF+J),J=1,LN),(AREF(IREF,K),K=K1,K2)
      CALL XFLUSH(6)
          END IF
50      CONTINUE
        WRITE(6,*)
      CALL XFLUSH(6)
100   CONTINUE
      WRITE(6,*)
      CALL XFLUSH(6)
      CALL QEXIT('REFCI')
      RETURN
      END
