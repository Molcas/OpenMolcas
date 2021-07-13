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
      SUBROUTINE STEPVEC(ICL,IOP,NCL,NOP,ISPIN,NORB,IWALK)
!
!     PURPOSE: A SPIN-COUPLED CSF IS SPECIFIED BY NCL CLOSED SHELL
!              AND NOP OPEN SHELL AND OCCUPATION VECTORS ICL AND IOP,
!              RESPECTIVELY. THE SPIN COUPLING IS STORED IN THE
!              VECTOR ISPIN. TRANSLATE THESE DATA INTO THE
!              CORRESPONDING STEP VECTOR.
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
      DIMENSION ICL(*),IOP(*),ISPIN(*)
      DIMENSION IWALK(*)
      Logical Test
!
      NXTCL = 1
      NXTOP = 1
      DO IORB = 1,NORB

        Test=NXTOP.LE.NOP
        If (Test) Test=IORB.EQ.IOP(NXTOP)

        IF(NXTCL.LE.NCL.AND.IORB.EQ.ICL(NXTCL) ) THEN

          IWALK(IORB) = 3
          NXTCL =NXTCL + 1

        ELSE IF( Test ) THEN

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

      END DO
!
!     EXIT
!
      RETURN
      END
