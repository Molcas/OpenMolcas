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
      SUBROUTINE STEPVEC(ICL,IOP,NCL,NOP,ISPIN,NORB,IWALK)
C
C     PURPOSE: A SPIN-COUPLED CSF IS SPECIFIED BY NCL CLOSED SHELL
C              AND NOP OPEN SHELL AND OCCUPATION VECTORS ICL AND IOP,
C              RESPECTIVELY. THE SPIN COUPLING IS STORED IN THE
C              VECTOR ISPIN. TRANSLATE THESE DATA INTO THE
C              CORRESPONDING STEP VECTOR.
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION ICL(*),IOP(*),ISPIN(*)
      DIMENSION IWALK(*)
C
      NXTCL = 1
      NXTOP = 1
      DO 100 IORB = 1,NORB
        IF(NXTCL.LE.NCL.AND.IORB.EQ.ICL(NXTCL) ) THEN
          IWALK(IORB) = 3
          NXTCL =NXTCL + 1
        ELSE IF(NXTOP.LE.NOP.AND.IORB.EQ.IOP(NXTOP) ) THEN
          IF(ISPIN(NXTOP).EQ.1) THEN
            IDELSP = 1
          ELSE
            IDELSP = -1
          END IF
          IF(IDELSP.EQ.1) THEN
             IWALK(IORB) = 1
          ELSE IF(IDELSP.EQ.-1) THEN
             IWALK(IORB) = 2
          END IF
          NXTOP = NXTOP + 1
        ELSE
          IWALK(IORB) = 0
        END IF
100   CONTINUE
C
C     EXIT
C
      RETURN
      END
