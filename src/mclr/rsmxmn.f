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
      SUBROUTINE RSMXMN(MAXEL,MINEL,NORB1,NORB2,NORB3,NEL,
     *                  MIN1,MAX1,MIN3,MAX3,IPRNT)
*
* Construct accumulated MAX and MIN arrays for a RAS set of strings
*
      IMPLICIT REAL*8 ( A-H,O-Z)
      DIMENSION  MINEL(*),MAXEL(*)
*
*. accumulated max and min in each of the three spaces
*. ( required max and min at final orbital in each space )
      MIN1A = MIN1
      MAX1A = MAX1
*
      MIN2A = NEL - MAX3
      MAX2A = NEL - MIN3
*
      MIN3A = NEL
      MAX3A = NEL
*
      NORB = NORB1 + NORB2 + NORB3
      DO 100 IORB = 1, NORB
        IF(IORB .LE. NORB1 ) THEN
          MINEL(IORB) = MAX(MIN1A+IORB-NORB1,0)
          MAXEL(IORB) = MIN(IORB,MAX1A)
        ELSE IF ( NORB1.LT.IORB .AND. IORB.LE.(NORB1+NORB2)) THEN
          MINEL(IORB) = MAX(MIN2A+IORB-NORB1-NORB2,0)
          IF(NORB1 .GT. 0 )
     *    MINEL(IORB) = MAX(MINEL(IORB),MINEL(NORB1))
          MAXEL(IORB) = MIN(IORB,MAX2A)
        ELSE IF ( IORB .GT. NORB1 + NORB2 ) THEN
          MINEL(IORB) = MAX(MIN3A+IORB-NORB,0)
          IF(NORB1+NORB2 .GT. 0 )
     *    MINEL(IORB) = MAX(MINEL(IORB),MINEL(NORB1+NORB2))
          MAXEL(IORB) = MIN(IORB,MAX3A)
        END IF
100   CONTINUE
*
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(IPRNT)
      END
