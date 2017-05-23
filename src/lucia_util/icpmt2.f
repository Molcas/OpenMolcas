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
      SUBROUTINE ICPMT2(     AIN,    AOUT,    NINR,    NINC,   NOUTR,
     &                     NOUTC,   IZERO)
*
* Copy INTEGER matrix AIN to AOUT . Dimensions can differ
*
* If IZERO .ne. 0 , AOUT is zeroed  first
      IMPLICIT REAL*8           (A-H,O-Z)
*. Input
      INTEGER AIN(NINR,NINC)
*. Output
      INTEGER AOUT(NOUTR,NOUTC)
*
      IF(IZERO.NE.0) CALL ISETVC(AOUT,0,NOUTR*NOUTC)
      DO 100 J = 1, NINC
       CALL ICOPVE(AIN(1,J),AOUT(1,J),NINR)
  100 CONTINUE
*
      RETURN
      END
