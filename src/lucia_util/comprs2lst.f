************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1996, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE COMPRS2LST(     I1,    XI1,     N1,     I2,    XI2,
     &                           N2,   NKIN,  NKOUT)
*
* Two lists of excitations/annihilations/creations are given.
* COmpress to common nonvanishing entries
*
* Jeppe Olsen, November 1996
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION I1(NKIN,N1),XI1(NKIN,N1)
      DIMENSION I2(NKIN,N2),XI2(NKIN,N2)
*
      NKOUT = 0
      DO K = 1, NKIN
        I1ACT  = 0
        DO I = 1, N1
          IF(I1(K,I).NE.0) I1ACT = 1
        END DO
        I2ACT = 0
        DO I = 1, N2
          IF(I2(K,I).NE.0) I2ACT = 1
        END DO
        IF(I1ACT.EQ.1.AND.I2ACT.EQ.1) THEN
          NKOUT = NKOUT + 1
          IF(NKOUT.NE.K) THEN
            DO I = 1, N1
               I1(NKOUT,I) = I1(K,I)
              XI1(NKOUT,I) =XI1(K,I)
            END DO
            DO I = 1, N2
               I2(NKOUT,I) = I2(K,I)
              XI2(NKOUT,I) =XI2(K,I)
            END DO
          END IF
        END IF
      END DO
*
      RETURN
      END
