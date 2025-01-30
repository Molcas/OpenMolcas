!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      SUBROUTINE MSSTRN_LUCIA(INSTRN,UTSTRN,NOPEN,IPRCSF)
!
! A STRING IS GIVEN IN FORM A SEQUENCE OF ZEROES
! AND ONE ' S
!
! REINTERPRET THIS AS :
!
! 1 : THE INPUT STRING IS A DETERMINANT AND THE
!     1'S INDICATE ALPHA ELECTRONS AND THE
!     0'S INDICATE BETA ELECTRONS .
!     UTSTRN IS THE MS-VALUES ATE EACH VERTEX
!
! 2 : THE INPUT STRING IS A CSF GIVEN IN A
!     BRANCHING DIAGRAM, WHERE
!     1'S INDICATE UPWARDS SPIN COUPLEING
!     WHILE THE 0'S INDICATES DOWNWARDS SPIN COUPLING ,
!     REEXPRESS THIS AS S VALUES OF ALL COUPLINGS
!
! THE TWO PROCEDURES ARE IDENTICAL .

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION INSTRN(NOPEN),UTSTRN(NOPEN)
!
      UTSTRN(1) = DBLE(INSTRN(1)) - 0.5D0
      DO 10 IOPEN = 2, NOPEN
        UTSTRN(IOPEN) = UTSTRN(IOPEN-1) +DBLE(INSTRN(IOPEN))-0.5D0
   10 CONTINUE
!
      NTEST = 0
      NTEST = MAX(NTEST,IPRCSF)
      IF(NTEST.GE.10) THEN
         WRITE(6,*) ' ... Output from MSSTRN '
         WRITE(6,*) ' INSTRN AND UTSTRN'
         CALL IWRTMA(INSTRN,1,NOPEN,1,NOPEN)
         CALL WRTMAT(UTSTRN,1,NOPEN,1,NOPEN)
       END IF
!
      RETURN
      END
