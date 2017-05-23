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
      SUBROUTINE NXTBLK_MCLR(IATP,IBTP,IASM,NOCTPA,NOCTPB,NSMST,IBLTP,
     &                  IDC,NONEW,IOCOC)
*
* Obtain allowed block following IATP IBTP IASM
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER IBLTP(*)
      INTEGER IOCOC(NOCTPA,NOCTPB)
*
C     CALL QENTER('NXTBLK')
*.Initialize
      ISM = IASM
      IA = IATP
      IB = IBTP
      NONEW = 0
*.Next block
 1000 CONTINUE
      IF(IB.LT.NOCTPB) THEN
        IB = IB + 1
      ELSE
        IB = 1
        IF(IA.LT.NOCTPA) THEN
          IA = IA+ 1
        ELSE
          IA = 1
          IF(ISM.LT.NSMST) THEN
            ISM = ISM + 1
          ELSE
            NONEW = 1
          END IF
        END IF
      END IF
      IF(NONEW.EQ.1) GOTO 1001
*. Should this block be included
      IF(IDC.NE.1.AND.IBLTP(ISM).EQ.0) GOTO 1000
      IF(IDC.NE.1.AND.IBLTP(ISM).EQ.2.AND.IA.LT.IB) GOTO 1000
      IF(IOCOC(IA,IB).EQ.0) GOTO 1000
 1001 CONTINUE
*
      IATP = IA
      IBTP = IB
      IASM = ISM
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,'(A,4I4)')
     &  ' NXTBLK : ISM IA IB NONEW ', IASM,IA,IB,NONEW
      END IF
*
C     CALL QEXIT('NXTBLK')
      RETURN
      END
