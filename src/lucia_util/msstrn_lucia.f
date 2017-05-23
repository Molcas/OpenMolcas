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
      SUBROUTINE MSSTRN_LUCIA(INSTRN,UTSTRN,NOPEN,IPRCSF)
C
C A STRING IS GIVEN IN FORM A SEQUENCE OF ZEROES
C AND ONE ' S
C
C REINTERPRET THIS AS :
C
C 1 : THE INPUT STRING IS A DETERMINANT AND THE
C     1'S INDICATE ALPHA ELECTRONS AND THE
C     0'S INDICATE BETA ELECTRONS .
C     UTSTRN IS THE MS-VALUES ATE EACH VERTEX
C
C 2 : THE INPUT STRING IS A CSF GIVEN IN A
C     BRANCHING DIAGRAM, WHERE
C     1'S INDICATE UPWARDS SPIN COUPLEING
C     WHILE THE 0'S INDICATES DOWNWARDS SPIN COUPLING ,
C     REEXPRESS THIS AS S VALUES OF ALL COUPLINGS
C
C THE TWO PROCEDURES ARE IDENTICAL .

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION INSTRN(NOPEN),UTSTRN(NOPEN)
C
      UTSTRN(1) = DBLE(INSTRN(1)) - 0.5D0
      DO 10 IOPEN = 2, NOPEN
        UTSTRN(IOPEN) = UTSTRN(IOPEN-1) +DBLE(INSTRN(IOPEN))-0.5D0
   10 CONTINUE
C
      NTEST = 0
      NTEST = MAX(NTEST,IPRCSF)
      IF(NTEST.GE.10) THEN
         WRITE(6,*) ' ... Output from MSSTRN '
         WRITE(6,*) ' INSTRN AND UTSTRN'
         CALL IWRTMA(INSTRN,1,NOPEN,1,NOPEN)
         CALL WRTMAT(UTSTRN,1,NOPEN,1,NOPEN)
       END IF
C
      RETURN
      END
