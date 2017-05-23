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
      SUBROUTINE CSFDET_LUCIA(  NOPEN,   IDET,   NDET,   ICSF,   NCSF,
     &                            CDC,   WORK, PSSIGN, IPRCSF)
*
* Expand csf's in terms of combinations with
* the use of the Graebenstetter method ( I.J.Q.C.10,P142(1976) )
*
* Input :
*         NOPEN : NUMBER OF OPEN ORBITALS
*         IDET  : OCCUPATION OF combinations
*         NDET  : NUMBER OF combinations
*         ICSF  : INTERMEDIATE SPIN COUPLINGS OF
*                 CSF'S IN BRANCHING DIAGRAM
* Output :
*         CDC :  NDET X NCSF MATRIX
*                GIVING EXPANSION FROM COMB'S TO CSF,S
*                CSF BASIS = Comb basis *CDC
* Scratch :
*          WORK ,SHOULD AT LEAST BE ???
*
* If combinations are use ( signaled by PSSIGN .ne. 0 )
* the factors are multiplies with sqrt(2), corresponding to
* a combination being 1/sqrt(2) times the sum or difference of two
* determinants
*
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IDET(NOPEN,NDET),ICSF(NOPEN,NCSF)
      DIMENSION CDC(NDET,NCSF)
      DIMENSION WORK(*)

      NTEST = 0
      NTEST = MAX(IPRCSF,NTEST)
      IF(PSSIGN.EQ.0.0D0) THEN
       CMBFAC = 1.0D0
      ELSE
       CMBFAC = SQRT(2.0D0)
      END IF
C
      KLFREE = 1
      KLMDET = KLFREE
      KLFREE = KLMDET + NDET * NOPEN
      KLSCSF =  KLFREE
      KLFREE = KLSCSF + NOPEN
C
C.. OBTAIN INTERMEDIATE VALUES OF MS FOR ALL DETERMINANTS
      DO 10 JDET = 1, NDET
        CALL MSSTRN_LUCIA(IDET(1,JDET),WORK(KLMDET+(JDET-1)*NOPEN),
     &        NOPEN,IPRCSF)
   10 CONTINUE
C
      DO 1000 JCSF = 1, NCSF
       IF( NTEST .GE. 105 ) WRITE(6,*) ' ....Output for CSF ',JCSF
C
C OBTAIN INTERMEDIATE COUPLINGS FOR CSF
      CALL MSSTRN_LUCIA(ICSF(1,JCSF),WORK(KLSCSF),NOPEN,IPRCSF)
C
      DO 900 JDET = 1, NDET
C EXPANSION COEFFICIENT OF DETERMINANT JDET FOR CSF JCSF
      COEF = 1.0D0
      SIGN = 1.0D0
      JDADD = (JDET-1)*NOPEN
      DO 700 IOPEN = 1, NOPEN
C
C + + CASE
        IF(ICSF(IOPEN,JCSF).EQ.1.AND.IDET(IOPEN,JDET).EQ.1) THEN
          COEF = COEF *
     &    (WORK(KLSCSF-1+IOPEN)+WORK(KLMDET-1+JDADD+IOPEN) )
     &    /(2.0D0*WORK(KLSCSF-1+IOPEN) )
        ELSE IF(ICSF(IOPEN,JCSF).EQ.1.AND.IDET(IOPEN,JDET).EQ.0) THEN
C + - CASE
          COEF = COEF *
     &    (WORK(KLSCSF-1+IOPEN)-WORK(KLMDET-1+JDADD+IOPEN) )
     &    /(2.0D0*WORK(KLSCSF-1+IOPEN) )
        ELSE IF(ICSF(IOPEN,JCSF).EQ.0.AND.IDET(IOPEN,JDET).EQ.1) THEN
C - + CASE
          COEF = COEF *
     &    (WORK(KLSCSF-1+IOPEN)-WORK(KLMDET-1+JDADD+IOPEN) +1.0D0)
     &    /(2.0D0*WORK(KLSCSF-1+IOPEN)+2.0D0 )
          SIGN  = - SIGN
        ELSE IF(ICSF(IOPEN,JCSF).EQ.0.AND.IDET(IOPEN,JDET).EQ.0) THEN
C - - CASE
          COEF = COEF *
     &    (WORK(KLSCSF-1+IOPEN)+WORK(KLMDET-1+JDADD+IOPEN) +1.0D0)
     &    /(2.0D0*WORK(KLSCSF-1+IOPEN)+2.0D0 )
        END IF
  700 CONTINUE
       CDC(JDET,JCSF) = SIGN * CMBFAC * SQRT(COEF)
  900 CONTINUE
 1000 CONTINUE
C
      IF( NTEST .GE. 5) THEN
        WRITE(6,*)
        WRITE(6,'(A,2I2)')
     &  '  The CDC array for  NOPEN ',NOPEN
        WRITE(6,*) ' NDET, NCSF = ', NDET,NCSF
        WRITE(6,*)
        CALL WRTMAT(CDC,NDET,NCSF,NDET,NCSF)
       END IF
C
      RETURN
      END
