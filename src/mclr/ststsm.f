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

      SUBROUTINE STSTSM_MCLR(STSTSX,STSTDX,NSMST)
*
* construct  STSTSX and STSTDX giving
* symmetry of sx (dx) connecting two given string symmetries
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER STSTSX(NSMST,NSMST),STSTDX(NSMST,NSMST)
*
      DO 100 ILSTSM = 1, NSMST
        DO 50 IRSTSM = 1, NSMST
          CALL SYMCOM_MCLR(1,5,ISXSM,IRSTSM,ILSTSM)
          CALL SYMCOM_MCLR(1,6,IDXSM,IRSTSM,ILSTSM)
          STSTSX(ILSTSM,IRSTSM) = ISXSM
          STSTDX(ILSTSM,IRSTSM) = IDXSM
   50   CONTINUE
  100 CONTINUE
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' STSTSM : STSTSX, STSTDX '
        CALL IWRTMA(STSTSX,NSMST,NSMST,NSMST,NSMST)
        CALL IWRTMA(STSTDX,NSMST,NSMST,NSMST,NSMST)
      END IF
*
      RETURN
      END
